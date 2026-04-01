"""
UHA Zero-Bias Test: Gaia/SH0ES vs Planck ξ comparison
=======================================================
Author: Eric D. Martin / All Your Baseline LLC
Date: 2026-03-24

Test the core UHA claim: ξ = d_c / d_H is invariant across H₀ values.

For objects where d_c is derived from REDSHIFT (not geometric distance):
  d_c(z, H₀) = (c / H₀) × integral_0^z dz' / E(z', Ω)
  d_H(H₀)    = c / H₀
  ξ(z, H₀)   = d_c / d_H = integral_0^z dz' / E(z', Ω)   [H₀ cancels exactly]

The RESIDUAL |Δξ| between SH0ES and Planck cosmologies comes ONLY from the
different Ω_m values (0.334 vs 0.315), not from H₀.

Data: Pantheon+SH0ES.dat — 1,701 SNe Ia with zCMB, MU_SH0ES, IS_CALIBRATOR

This script measures the geometric (Ωm-driven) residual in ξ-space between
the SH0ES-like and Planck-like frameworks. The Ωm difference (0.334 vs 0.315)
is the sole source of any |Δξ| residual; H₀ cancels exactly in d_c/d_H.
"""

import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import matplotlib.pyplot as plt
import os

# ---------------------------------------------------------------------------
# Cosmological parameters
# ---------------------------------------------------------------------------
# SH0ES 2022 (Riess et al. 2022)
SH0ES = dict(H0=73.04, Om0=0.334)

# Planck 2018 (Planck Collaboration 2020)
PLANCK = dict(H0=67.4, Om0=0.315)

# Speed of light
c_kms = 299792.458  # km/s

# Horizon radius d_H = c / H₀
d_H_SH0ES  = c_kms / SH0ES['H0']   # Mpc
d_H_PLANCK = c_kms / PLANCK['H0']  # Mpc

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
data_path = os.path.join(os.path.dirname(__file__), '..', 'cosmological-data',
                         'Pantheon+_Data', '4_DISTANCES_AND_COVAR', 'Pantheon+SH0ES.dat')
df = pd.read_csv(data_path, sep=r'\s+', comment='#')

print(f"Loaded {len(df)} rows from Pantheon+SH0ES.dat")
print(f"Calibrators (IS_CALIBRATOR=1): {(df['IS_CALIBRATOR'] == 1).sum()}")
print(f"Cosmological SNe:              {(df['IS_CALIBRATOR'] == 0).sum()}")

# ---------------------------------------------------------------------------
# Build astropy cosmology objects
# ---------------------------------------------------------------------------
cosmo_shoes  = FlatLambdaCDM(H0=SH0ES['H0'],  Om0=SH0ES['Om0'])
cosmo_planck = FlatLambdaCDM(H0=PLANCK['H0'], Om0=PLANCK['Om0'])

# ---------------------------------------------------------------------------
# Compute ξ from redshift for each object
# Method: d_c(z) from each cosmology, then ξ = d_c / d_H
# H₀ cancels exactly in the ratio — residual is from Ω_m only
# ---------------------------------------------------------------------------
z = df['zCMB'].values
z_ok = z > 0.005  # exclude very nearby objects where peculiar velocities dominate

df_cos = df[z_ok].copy()
z_cos = df_cos['zCMB'].values

print(f"\nUsing {len(df_cos)} objects with z > 0.005")

# Comoving distances [Mpc]
d_c_shoes  = cosmo_shoes.comoving_distance(z_cos).to(u.Mpc).value
d_c_planck = cosmo_planck.comoving_distance(z_cos).to(u.Mpc).value

# ξ = d_c / d_H (dimensionless, horizon-normalized)
xi_shoes  = d_c_shoes  / d_H_SH0ES
xi_planck = d_c_planck / d_H_PLANCK

# Key invariance check: since d_c ∝ 1/H₀ and d_H ∝ 1/H₀, both cancel.
# The ratio ξ depends only on (z, Ω_m) — H₀-free.
# Verify: xi_shoes_from_function ≈ f(z, Ω_m=0.334), xi_planck ≈ f(z, Ω_m=0.315)
#
# For z << 1: xi ≈ z (independent of all cosmological parameters)
# Residual Δξ comes from nonlinear Ω_m terms at higher z

delta_xi = xi_shoes - xi_planck
abs_delta = np.abs(delta_xi)

# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------
print("\n" + "="*60)
print("UHA Zero-Bias Test Results")
print("="*60)
print(f"\nH₀ values:  SH0ES = {SH0ES['H0']} km/s/Mpc,  Planck = {PLANCK['H0']} km/s/Mpc")
print(f"Ω_m values: SH0ES = {SH0ES['Om0']},         Planck = {PLANCK['Om0']}")
print(f"H₀ tension: {(SH0ES['H0']-PLANCK['H0'])/PLANCK['H0']*100:.1f}% ({(SH0ES['H0']-PLANCK['H0']):.2f} km/s/Mpc)")
print(f"\nξ statistics (SH0ES cosmology):")
print(f"  Mean ξ:   {xi_shoes.mean():.6f}")
print(f"  Range:    [{xi_shoes.min():.6f}, {xi_shoes.max():.6f}]")
print(f"\n|Δξ| = |ξ_SH0ES - ξ_Planck|:")
print(f"  Mean:     {abs_delta.mean():.2e}")
print(f"  Median:   {np.median(abs_delta):.2e}")
print(f"  Max:      {abs_delta.max():.2e}")
print(f"  Std:      {abs_delta.std():.2e}")
print(f"\nFraction with |Δξ| < 1e-2: {(abs_delta < 1e-2).mean()*100:.1f}%")
print(f"Fraction with |Δξ| < 1e-3: {(abs_delta < 1e-3).mean()*100:.1f}%")

# ξ-space residual as fraction of H₀ gap
raw_tension = abs(SH0ES['H0'] - PLANCK['H0']) / PLANCK['H0']
xi_tension  = abs_delta.mean() / np.mean((xi_shoes + xi_planck) / 2)
xi_fraction_of_gap = xi_tension / raw_tension * 100
print(f"\nH₀ fractional gap (SH0ES–Planck):     {raw_tension*100:.1f}%")
print(f"ξ fractional residual (ratio-of-means): {xi_tension*100:.4f}%")
print(f"ξ residual as fraction of H₀ gap:      {xi_fraction_of_gap:.1f}%")

# ---------------------------------------------------------------------------
# Low-z limit: ξ ≈ z (pure invariance regime)
# ---------------------------------------------------------------------------
low_z = z_cos < 0.1
print(f"\nLow-z objects (z < 0.1): {low_z.sum()}")
print(f"  Mean |Δξ| at low z:    {abs_delta[low_z].mean():.2e}  (expected ≈ 0 as z→0)")

# High-z limit: Ω_m difference matters more
high_z = z_cos > 0.5
print(f"\nHigh-z objects (z > 0.5): {high_z.sum()}")
print(f"  Mean |Δξ| at high z:    {abs_delta[high_z].mean():.2e}")

# ---------------------------------------------------------------------------
# Comparison: ξ vs raw distance modulus scatter
# ---------------------------------------------------------------------------
# In ξ space: the H₀ factor in distance ladder is removed
# The remaining scatter is the physical Ω_m residual
print(f"\n" + "="*60)
print("Interpretation")
print("="*60)
print(f"""
The nominal H₀ fractional gap (SH0ES–Planck) is {raw_tension*100:.1f}%.
The ξ-space residual attributable to the Ωm difference (0.334 vs 0.315)
is {xi_tension*100:.3f}% (ratio-of-means: ⟨|Δξ|⟩ / ⟨ξ_mid⟩).

This represents {xi_fraction_of_gap:.1f}% of the nominal H₀ gap.
The remainder reflects other physical or systematic sources and is
not quantified by this diagnostic.

ξ-space diagnostic:
  ⟨|Δξ|⟩     = {abs_delta.mean():.2e}
  ⟨ξ_mid⟩    = {np.mean((xi_shoes+xi_planck)/2):.4f}
  ratio       = {xi_tension*100:.4f}%
  H₀ gap      = {raw_tension*100:.1f}%
  fraction    = {xi_fraction_of_gap:.1f}%
""")

# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------
out_dir = os.path.join(os.path.dirname(__file__), '..', 'papers')
results = pd.DataFrame({
    'z': z_cos,
    'xi_shoes': xi_shoes,
    'xi_planck': xi_planck,
    'delta_xi': delta_xi,
    'abs_delta_xi': abs_delta,
})
results_path = os.path.join(out_dir, 'gaia_zero_bias_results.csv')
results.to_csv(results_path, index=False)
print(f"Results saved to {results_path}")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

axes[0].scatter(z_cos, xi_shoes,  s=1, alpha=0.3, label=f'SH0ES H₀={SH0ES["H0"]}', c='steelblue')
axes[0].scatter(z_cos, xi_planck, s=1, alpha=0.3, label=f'Planck H₀={PLANCK["H0"]}', c='darkorange')
axes[0].set_xlabel('Redshift z'); axes[0].set_ylabel('ξ = d_c / d_H')
axes[0].set_title('ξ vs z: SH0ES vs Planck')
axes[0].legend(markerscale=5)

axes[1].scatter(z_cos, abs_delta, s=1, alpha=0.3, c='crimson')
axes[1].set_xlabel('Redshift z'); axes[1].set_ylabel('|Δξ|')
axes[1].set_title('|Δξ| across cosmologies')
axes[1].set_yscale('log')

axes[2].hist(np.log10(abs_delta[abs_delta > 0]), bins=50, color='purple', alpha=0.7)
axes[2].set_xlabel('log₁₀|Δξ|'); axes[2].set_ylabel('Count')
axes[2].set_title(f'Distribution of |Δξ| (n={len(abs_delta)})')
axes[2].axvline(np.log10(raw_tension), color='red', linestyle='--',
                label=f'Raw H₀ tension ({raw_tension*100:.1f}%)')
axes[2].legend()

plt.tight_layout()
plot_path = os.path.join(out_dir, 'gaia_zero_bias_plot.png')
plt.savefig(plot_path, dpi=150, bbox_inches='tight')
print(f"Plot saved to {plot_path}")
plt.close()
