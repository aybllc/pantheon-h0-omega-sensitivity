"""
Bootstrap V1 — ΔH₀ stability under Cepheid-anchored refit.
1000 resamples of Hubble-flow SNe. Reports mean, std, S/N.
"""
import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from scipy.optimize import minimize_scalar
import os

np.random.seed(42)

data_path = os.path.join(os.path.dirname(__file__), '..', 'Pantheon+_Data',
                         '4_DISTANCES_AND_COVAR', 'Pantheon+SH0ES.dat')
df = pd.read_csv(data_path, sep=r'\s+', comment='#')

cal = df[df['IS_CALIBRATOR'] == 1].copy()
hf  = df[(df['IS_CALIBRATOR'] == 0) & (df['zCMB'] > 0.005)].copy()

z_cal     = cal['zCMB'].values
mu_cal    = cal['CEPH_DIST'].values
sigma_cal = cal['MU_SH0ES_ERR_DIAG'].values

z_hf     = hf['zCMB'].values
mu_hf    = hf['MU_SH0ES'].values
sigma_hf = hf['MU_SH0ES_ERR_DIAG'].values

H0_ref = 73.04
Om_ref = 0.334

def mu_theory(z, H0, Om0):
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    return 5.0 * np.log10(cosmo.luminosity_distance(z).to(u.Mpc).value) + 25.0

# Anchor M at reference cosmology from calibrators
w_cal   = 1.0 / sigma_cal**2
mu_th_cal_ref = mu_theory(z_cal, H0_ref, Om_ref)
M_fixed = np.sum(w_cal * (mu_cal - mu_th_cal_ref)) / np.sum(w_cal)

def chi2_hf(H0, Om0, z, mu, sig):
    resid = mu - (mu_theory(z, H0, Om0) + M_fixed)
    return np.sum((resid / sig)**2)

def best_H0(Om0, z, mu, sig, bracket=(68.0, 82.0)):
    H0_scan  = np.linspace(bracket[0], bracket[1], 71)
    chi2_scan = np.array([chi2_hf(h, Om0, z, mu, sig) for h in H0_scan])
    H0_init  = H0_scan[np.argmin(chi2_scan)]
    lo = max(bracket[0], H0_init - 2.0)
    hi = min(bracket[1], H0_init + 2.0)
    res = minimize_scalar(chi2_hf, bounds=(lo, hi), method='bounded',
                          args=(Om0, z, mu, sig), options={'xatol': 1e-4})
    return res.x

Om_shoes  = 0.334
Om_planck = 0.315
N_boot    = 1000
n_hf      = len(hf)

print(f"Bootstrap V1: {N_boot} resamples, n_hf={n_hf}")
print(f"M_fixed = {M_fixed:.5f}")

delta_H0_boot = np.zeros(N_boot)

for i in range(N_boot):
    idx = np.random.choice(n_hf, size=n_hf, replace=True)
    zb  = z_hf[idx]
    mub = mu_hf[idx]
    sb  = sigma_hf[idx]

    h_s = best_H0(Om_shoes,  zb, mub, sb)
    h_p = best_H0(Om_planck, zb, mub, sb)
    delta_H0_boot[i] = h_s - h_p

    if (i+1) % 100 == 0:
        print(f"  {i+1}/{N_boot}  running mean={delta_H0_boot[:i+1].mean():.4f}  "
              f"std={delta_H0_boot[:i+1].std():.4f}")

out_path = os.path.join(os.path.dirname(__file__), 'bootstrap_deltaH0.npy')
np.save(out_path, delta_H0_boot)

mean_dH0 = delta_H0_boot.mean()
std_dH0  = delta_H0_boot.std()
sn       = mean_dH0 / std_dH0
pct_zero = np.mean(delta_H0_boot <= 0) * 100

print("\n" + "="*50)
print("BOOTSTRAP V1 RESULT")
print("="*50)
print(f"  Mean ΔH₀   = {mean_dH0:+.4f} km/s/Mpc")
print(f"  Std  ΔH₀   = {std_dH0:.4f} km/s/Mpc")
print(f"  S/N        = {sn:.2f}")
print(f"  P(ΔH₀≤0)  = {pct_zero:.1f}%")
print(f"  Saved → {out_path}")

if sn >= 1.0 and pct_zero < 16.0:
    print("\nPASS: ΔH₀ is stable and consistently positive.")
else:
    print("\nFAIL: ΔH₀ not stable — review before claiming.")
