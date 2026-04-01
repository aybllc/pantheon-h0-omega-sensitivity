"""
Sensitivity of Pantheon+ Hubble-flow H₀ to assumed matter density Ωm.
======================================================================
Author: Eric D. Martin / All Your Baseline LLC
Date: 2026-04-01

Correct pipeline (v2):
  Anchors the absolute magnitude M using Cepheid calibrators (IS_CALIBRATOR=1,
  CEPH_DIST) to break the H₀-M degeneracy, then fits H₀ from the Hubble-flow
  sample at fixed Ωm.

Step 1 — Compute M from calibrators at reference (H₀=73.04, Ωm=0.334)
Step 2 — Fix M
Step 3 — Fit H₀ from Hubble-flow SNe at Ωm=0.334  [recover ~73.04]
Step 4 — Fit H₀ from Hubble-flow SNe at Ωm=0.315
Step 5 — Report ΔH₀ and H₀(Ωm) scan
"""

import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
import os

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
data_path = os.path.join(os.path.dirname(__file__), '..', 'Pantheon+_Data',
                         '4_DISTANCES_AND_COVAR', 'Pantheon+SH0ES.dat')
df = pd.read_csv(data_path, sep=r'\s+', comment='#')

# Calibrators: Cepheid-anchored hosts, z << 1
cal  = df[df['IS_CALIBRATOR'] == 1].copy()

# Hubble-flow sample: non-calibrators, z > 0.005
hf_mask = (df['IS_CALIBRATOR'] == 0) & (df['zCMB'] > 0.005)
hf   = df[hf_mask].copy()

z_hf      = hf['zCMB'].values
mu_hf     = hf['MU_SH0ES'].values
sigma_hf  = hf['MU_SH0ES_ERR_DIAG'].values

z_cal     = cal['zCMB'].values
mu_cal    = cal['CEPH_DIST'].values        # geometric, H₀/Ωm independent
sigma_cal = cal['MU_SH0ES_ERR_DIAG'].values

print(f"Calibrators:       {len(cal)} SNe  (z range: {z_cal.min():.4f} – {z_cal.max():.4f})")
print(f"Hubble-flow SNe:   {len(hf)}  (z range: {z_hf.min():.4f} – {z_hf.max():.3f})")

# ---------------------------------------------------------------------------
# Helper: theoretical distance modulus
# ---------------------------------------------------------------------------
def mu_theory(z, H0, Om0):
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    d_L   = cosmo.luminosity_distance(z).to(u.Mpc).value
    return 5.0 * np.log10(d_L) + 25.0

# ---------------------------------------------------------------------------
# Step 1 — Anchor M from calibrators
#
# For each calibrator: MU_obs = CEPH_DIST (geometric)
#                      Model:   μ_theory(z_cal, H0, Ωm) + M
# At the calibrators' low z (< 0.017), μ_theory barely depends on Ωm.
# We use the SH0ES reference cosmology to define M, then hold it fixed.
#
# M encodes: the offset between the Cepheid distance scale and the
# cosmological H₀ scale.  At fixed M, changing Ωm in the Hubble flow
# reveals how H₀ must shift to match the data.
# ---------------------------------------------------------------------------
H0_ref  = 73.04
Om_ref  = 0.334

mu_th_cal_ref = mu_theory(z_cal, H0_ref, Om_ref)

w_cal   = 1.0 / sigma_cal**2
M_fixed = np.sum(w_cal * (mu_cal - mu_th_cal_ref)) / np.sum(w_cal)

print(f"\nCalibrator-anchored M = {M_fixed:.5f} mag")
print(f"  (residual vs reference cosmology; should be near 0)")

# ---------------------------------------------------------------------------
# Step 2 — χ²(H₀) for Hubble-flow sample at fixed M and fixed Ωm
# ---------------------------------------------------------------------------
def chi2_hf(H0, Om0):
    """χ² at fixed M. H₀ is the only free parameter."""
    mu_th  = mu_theory(z_hf, H0, Om0) + M_fixed
    resid  = mu_hf - mu_th
    return np.sum((resid / sigma_hf)**2)

def best_H0(Om0, bracket=(65.0, 80.0)):
    """Robustly find H₀ minimizing χ²(H₀) at fixed Ωm."""
    # Grid scan first → identify approximate minimum
    H0_scan  = np.linspace(bracket[0], bracket[1], 151)
    chi2_scan = np.array([chi2_hf(h, Om0) for h in H0_scan])
    H0_init  = H0_scan[np.argmin(chi2_scan)]
    # Refine near the grid minimum
    lo = max(bracket[0], H0_init - 2.0)
    hi = min(bracket[1], H0_init + 2.0)
    res = minimize_scalar(chi2_hf, bounds=(lo, hi),
                          method='bounded', args=(Om0,),
                          options={'xatol': 1e-5})
    return res.x, res.fun

# ---------------------------------------------------------------------------
# Step 3/4 — Key measurements
# ---------------------------------------------------------------------------
Om_shoes  = 0.334
Om_planck = 0.315

print(f"\nFitting H₀ at Ωm = {Om_shoes} (SH0ES-like)...")
H0_shoes,  chi2_s = best_H0(Om_shoes)
print(f"  → H₀ = {H0_shoes:.3f} km/s/Mpc  (χ²/dof = {chi2_s/(len(hf)-1):.4f})")

print(f"\nFitting H₀ at Ωm = {Om_planck} (Planck-like)...")
H0_planck, chi2_p = best_H0(Om_planck)
print(f"  → H₀ = {H0_planck:.3f} km/s/Mpc  (χ²/dof = {chi2_p/(len(hf)-1):.4f})")

delta_H0     = H0_shoes - H0_planck
delta_H0_pct = delta_H0 / H0_planck * 100.0

print("\n" + "="*60)
print("RESULT 2: Direct H₀ Refit — Ωm Sensitivity")
print("="*60)
print(f"  Ωm = {Om_shoes}  →  H₀ = {H0_shoes:.3f} km/s/Mpc")
print(f"  Ωm = {Om_planck}  →  H₀ = {H0_planck:.3f} km/s/Mpc")
print(f"  ΔH₀ = {delta_H0:+.3f} km/s/Mpc  ({delta_H0_pct:+.2f}% of H₀)")

# ---------------------------------------------------------------------------
# Step 5 — H₀(Ωm) scan
# ---------------------------------------------------------------------------
print("\nRunning smooth H₀(Ωm) scan (17 points)...")
omega_grid = np.linspace(0.28, 0.36, 17)
H0_grid    = []
for Om0 in omega_grid:
    h, _ = best_H0(Om0)
    H0_grid.append(h)
    print(f"  Ωm = {Om0:.4f}  →  H₀ = {h:.3f}")
H0_grid = np.array(H0_grid)

# ---------------------------------------------------------------------------
# Cross-check with Result 1
# ---------------------------------------------------------------------------
xi_frac_pct = 0.57    # from gaia_zero_bias_test.py
H0_gap_pct  = abs((73.04 - 67.4) / 67.4 * 100)   # nominal SH0ES - Planck

print("\n" + "="*60)
print("CROSS-CHECK: ξ-space vs direct refit")
print("="*60)
print(f"  Nominal SH0ES – Planck H₀ gap:  {H0_gap_pct:.1f}%")
print(f"  ξ-space fractional residual:    {xi_frac_pct:.2f}%")
print(f"  Direct refit ΔH₀:               {abs(delta_H0_pct):.2f}%")

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out_dir  = os.path.dirname(__file__)
scan_df  = pd.DataFrame({'Om0': omega_grid, 'H0_best': H0_grid})
scan_path = os.path.join(out_dir, 'refit_H0_vs_Om.csv')
scan_df.to_csv(scan_path, index=False)
print(f"\nScan saved → {scan_path}")

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(omega_grid, H0_grid, 'k-o', lw=2, ms=5,
        label=r'Best-fit $H_0(\Omega_m)$ [Cepheid-anchored]')
ax.axvline(Om_shoes,  color='steelblue',  lw=1.5, ls='--',
           label=rf'SH0ES $\Omega_m={Om_shoes}$ → {H0_shoes:.2f}')
ax.axvline(Om_planck, color='darkorange', lw=1.5, ls='--',
           label=rf'Planck $\Omega_m={Om_planck}$ → {H0_planck:.2f}')
ax.axhspan(min(H0_shoes, H0_planck), max(H0_shoes, H0_planck),
           alpha=0.08, color='gray',
           label=rf'$\Delta H_0={delta_H0:+.2f}$ km/s/Mpc ({delta_H0_pct:+.1f}%)')
ax.set_xlabel(r'Assumed $\Omega_m$', fontsize=13)
ax.set_ylabel(r'Best-fit $H_0$ (km/s/Mpc)', fontsize=13)
ax.set_title(r'Sensitivity of Pantheon+ Hubble-flow $H_0$ to assumed $\Omega_m$'
             f'\n(Cepheid-anchored; n={len(hf)} Hubble-flow SNe)', fontsize=11)
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
plt.tight_layout()
fig_path = os.path.join(out_dir, 'refit_H0_vs_Om.png')
plt.savefig(fig_path, dpi=200, bbox_inches='tight')
print(f"Figure saved → {fig_path}")
plt.close()

print("\n" + "="*60)
print("PAPER RESULT SUMMARY (both diagnostics)")
print("="*60)
print(f"""
Result 1 — ξ-space (gaia_zero_bias_test.py):
  Mean fractional ξ residual = {xi_frac_pct:.2f}% across {len(hf)} Hubble-flow SNe Ia
  between SH0ES-like (Ωm=0.334) and Planck-like (Ωm=0.315) frameworks.

Result 2 — Direct Cepheid-anchored refit (this script):
  ΔH₀ = {delta_H0:+.3f} km/s/Mpc ({delta_H0_pct:+.2f}%)
  from Ωm = {Om_shoes} → {Om_planck}, with M fixed from {len(cal)} Cepheid calibrators.

Nominal H₀ gap: {H0_gap_pct:.1f}% ({73.04 - 67.4:.2f} km/s/Mpc)

Both diagnostics quantify framework dependence of Hubble-flow H₀.
""")
