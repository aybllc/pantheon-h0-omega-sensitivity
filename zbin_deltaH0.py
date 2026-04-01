"""
V2 — z-binned ΔH₀ analysis.
Split 1616 Hubble-flow SNe into 4 z-bins.
Compute H₀(Ωm=0.334) and H₀(Ωm=0.315) per bin.
Check uniformity of Ωm sensitivity across redshift.
"""
import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
import os

data_path = os.path.join(os.path.dirname(__file__), '..', 'Pantheon+_Data',
                         '4_DISTANCES_AND_COVAR', 'Pantheon+SH0ES.dat')
df = pd.read_csv(data_path, sep=r'\s+', comment='#')

cal = df[df['IS_CALIBRATOR'] == 1].copy()
hf  = df[(df['IS_CALIBRATOR'] == 0) & (df['zCMB'] > 0.005)].copy()

# Anchor M from calibrators at reference cosmology
H0_ref, Om_ref = 73.04, 0.334
def mu_theory(z, H0, Om0):
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    return 5.0 * np.log10(cosmo.luminosity_distance(z).to(u.Mpc).value) + 25.0

z_cal     = cal['zCMB'].values
mu_cal    = cal['CEPH_DIST'].values
sigma_cal = cal['MU_SH0ES_ERR_DIAG'].values
w_cal     = 1.0 / sigma_cal**2
M_fixed   = np.sum(w_cal * (mu_cal - mu_theory(z_cal, H0_ref, Om_ref))) / np.sum(w_cal)
print(f"M_fixed = {M_fixed:.5f}")

def chi2_hf(H0, Om0, z, mu, sig):
    resid = mu - (mu_theory(z, H0, Om0) + M_fixed)
    return np.sum((resid / sig)**2)

def best_H0(Om0, z, mu, sig, bracket=(68.0, 82.0)):
    H0_scan   = np.linspace(bracket[0], bracket[1], 71)
    chi2_scan = np.array([chi2_hf(h, Om0, z, mu, sig) for h in H0_scan])
    H0_init   = H0_scan[np.argmin(chi2_scan)]
    lo = max(bracket[0], H0_init - 2.0)
    hi = min(bracket[1], H0_init + 2.0)
    res = minimize_scalar(chi2_hf, bounds=(lo, hi), method='bounded',
                          args=(Om0, z, mu, sig), options={'xatol': 1e-4})
    return res.x

Om_shoes  = 0.334
Om_planck = 0.315

# Equal-count bins
z_hf     = hf['zCMB'].values
mu_hf    = hf['MU_SH0ES'].values
sig_hf   = hf['MU_SH0ES_ERR_DIAG'].values

edges = np.percentile(z_hf, [0, 25, 50, 75, 100])
print(f"\nz-bin edges (quartiles): {edges}")

rows = []
for i in range(4):
    mask = (z_hf >= edges[i]) & (z_hf < edges[i+1]) if i < 3 else (z_hf >= edges[i])
    zb, mub, sb = z_hf[mask], mu_hf[mask], sig_hf[mask]
    h_s = best_H0(Om_shoes,  zb, mub, sb)
    h_p = best_H0(Om_planck, zb, mub, sb)
    dH  = h_s - h_p
    rows.append({'bin': i+1, 'z_lo': edges[i], 'z_hi': edges[i+1],
                 'n': mask.sum(), 'H0_shoes': h_s, 'H0_planck': h_p, 'deltaH0': dH})
    print(f"  Bin {i+1}: z=[{edges[i]:.3f},{edges[i+1]:.3f}]  n={mask.sum():4d}  "
          f"H0(0.334)={h_s:.3f}  H0(0.315)={h_p:.3f}  ΔH0={dH:+.4f}")

res_df = pd.DataFrame(rows)
out_csv = os.path.join(os.path.dirname(__file__), 'zbin_deltaH0.csv')
res_df.to_csv(out_csv, index=False)

print("\n" + "="*55)
print("V2 z-BINNED RESULT")
print("="*55)
print(f"  ΔH₀ per bin: {res_df['deltaH0'].values}")
print(f"  Range:  [{res_df['deltaH0'].min():+.4f}, {res_df['deltaH0'].max():+.4f}]")
print(f"  Mean:   {res_df['deltaH0'].mean():+.4f}")
print(f"  Std:    {res_df['deltaH0'].std():.4f}")
all_same_sign = (res_df['deltaH0'] < 0).all() or (res_df['deltaH0'] > 0).all()
print(f"  Uniform sign: {all_same_sign}")
if all_same_sign:
    print("  PASS: ΔH₀ sign consistent across all z-bins.")
else:
    print("  NOTE: sign flip detected — sensitivity is z-dependent.")

# Figure
fig, ax = plt.subplots(figsize=(7, 4))
z_mids = [(r['z_lo'] + r['z_hi']) / 2 for _, r in res_df.iterrows()]
ax.bar(range(1, 5), res_df['deltaH0'], color='steelblue', alpha=0.8)
ax.axhline(0, color='k', lw=0.8, ls='--')
ax.axhline(res_df['deltaH0'].mean(), color='darkorange', lw=1.5, ls='--',
           label=f"mean={res_df['deltaH0'].mean():+.4f}")
ax.set_xticks(range(1, 5))
ax.set_xticklabels([f"z=[{r['z_lo']:.2f},{r['z_hi']:.2f}]\nn={r['n']}"
                    for _, r in res_df.iterrows()], fontsize=8)
ax.set_ylabel(r'$\Delta H_0$ (km/s/Mpc)', fontsize=12)
ax.set_title(r'$\Delta H_0 = H_0(\Omega_m=0.334) - H_0(\Omega_m=0.315)$ by z-bin', fontsize=11)
ax.legend(fontsize=9)
ax.grid(alpha=0.3, axis='y')
plt.tight_layout()
fig_path = os.path.join(os.path.dirname(__file__), 'zbin_deltaH0.png')
plt.savefig(fig_path, dpi=200, bbox_inches='tight')
print(f"\nSaved → {out_csv}")
print(f"Saved → {fig_path}")
plt.close()
