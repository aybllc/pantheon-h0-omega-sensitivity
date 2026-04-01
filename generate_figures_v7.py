"""
generate_figures_v7.py
Generate all new figures needed for hubble_om_sensitivity_v7.tex.

Figures produced (all saved to same directory as this script):
  bootstrap_histogram.png
  xi_vs_z_scatter.png
  xi_residual_cdf.png
  chi2_sanity_check.png
  chi2_comparison.png
  residual_vs_z.png
  calibrator_ceph_vs_mu.png
  z_distribution.png
  sigma_vs_z.png
  H0_scan_per_zbin.png
  H0_comparison_ladder.png
  ximid_vs_z_bins.png

Requires: bootstrap_deltaH0.npy, refit_H0_vs_Om.csv, zbin_deltaH0.csv
          (all in same directory), and Pantheon+SH0ES.dat (see data_path).
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from scipy.optimize import minimize_scalar

# ----------------------------------------------------------------
# Paths
# ----------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(HERE, '..', 'cosmological-data',
                         'Pantheon+_Data',
                         '4_DISTANCES_AND_COVAR',
                         'Pantheon+SH0ES.dat')

def fig_path(name):
    return os.path.join(HERE, name)

# ----------------------------------------------------------------
# Load data
# ----------------------------------------------------------------
print("Loading data ...")
df = pd.read_csv(DATA_PATH, sep=r'\s+', comment='#')
cal  = df[df['IS_CALIBRATOR'] == 1].copy()
hf   = df[(df['IS_CALIBRATOR'] == 0) & (df['zCMB'] > 0.005)].copy()
full = df[df['zCMB'] > 0.001].copy()   # all SNe for xi diagnostic

# Reference cosmology / anchor
H0_ref, Om_ref = 73.04, 0.334
Om_shoes, Om_planck = 0.334, 0.315

def mu_theory(z, H0, Om0):
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    return 5.0 * np.log10(cosmo.luminosity_distance(z).to(u.Mpc).value) + 25.0

def xi_func(z, Om0):
    cosmo = FlatLambdaCDM(H0=100.0, Om0=Om0)
    dc = cosmo.comoving_distance(z).to(u.Mpc).value
    dH = 2997.92458  # c/H0 at H0=100, in Mpc
    return dc / dH

# Cepheid-anchored M
z_cal     = cal['zCMB'].values
mu_cal    = cal['CEPH_DIST'].values
sig_cal   = cal['MU_SH0ES_ERR_DIAG'].values
w_cal     = 1.0 / sig_cal**2
M_fixed   = np.sum(w_cal * (mu_cal - mu_theory(z_cal, H0_ref, Om_ref))) / np.sum(w_cal)
print(f"  M_fixed = {M_fixed:.5f}")

z_hf  = hf['zCMB'].values
mu_hf = hf['MU_SH0ES'].values
sig_hf = hf['MU_SH0ES_ERR_DIAG'].values

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

H0_s = best_H0(Om_shoes,  z_hf, mu_hf, sig_hf)
H0_p = best_H0(Om_planck, z_hf, mu_hf, sig_hf)
print(f"  H0(0.334) = {H0_s:.3f}   H0(0.315) = {H0_p:.3f}   DH0 = {H0_s-H0_p:+.4f}")

# ----------------------------------------------------------------
# Load existing derived data
# ----------------------------------------------------------------
boot = np.load(os.path.join(HERE, 'bootstrap_deltaH0.npy'))
h0scan_df = pd.read_csv(os.path.join(HERE, 'refit_H0_vs_Om.csv'))
zbin_df   = pd.read_csv(os.path.join(HERE, 'zbin_deltaH0.csv'))

# xi values for all SNe
z_full = full['zCMB'].values
xi_s   = xi_func(z_full, Om_shoes)
xi_p   = xi_func(z_full, Om_planck)
xi_mid = (xi_s + xi_p) / 2.0
dxi    = np.abs(xi_s - xi_p)
frac   = dxi / xi_mid

STYLE = dict(dpi=200, bbox_inches='tight')
plt.rcParams.update({'font.size': 10})

# ----------------------------------------------------------------
# 1. Bootstrap histogram
# ----------------------------------------------------------------
print("Fig 1: bootstrap_histogram ...")
fig, ax = plt.subplots(figsize=(6, 4))
ax.hist(boot, bins=40, color='steelblue', alpha=0.85, edgecolor='white')
ax.axvline(boot.mean(), color='k', lw=1.5, ls='-', label=f'Mean = {boot.mean():.4f}')
ax.axvspan(np.percentile(boot, 16), np.percentile(boot, 84),
           alpha=0.15, color='steelblue', label='16th–84th pct')
ax.set_xlabel(r'$\Delta H_0$ (km s$^{-1}$ Mpc$^{-1}$)', fontsize=11)
ax.set_ylabel('Count', fontsize=11)
ax.set_title(r'Bootstrap distribution: $\Delta H_0 = H_0(0.334) - H_0(0.315)$', fontsize=11)
ax.legend(fontsize=9)
ax.grid(alpha=0.3, axis='y')
fig.tight_layout()
fig.savefig(fig_path('bootstrap_histogram.png'), **STYLE)
plt.close()

# ----------------------------------------------------------------
# 2. |Δξ/ξ| vs z scatter
# ----------------------------------------------------------------
print("Fig 2: xi_vs_z_scatter ...")
fig, ax = plt.subplots(figsize=(7, 4))
ax.scatter(z_full, frac * 100, s=1, c='grey', alpha=0.4, rasterized=True)
# running mean in log-z bins
log_edges = np.linspace(np.log10(z_full.min()), np.log10(z_full.max()), 25)
zmid_run, frac_run = [], []
for lo, hi in zip(log_edges[:-1], log_edges[1:]):
    mask = (np.log10(z_full) >= lo) & (np.log10(z_full) < hi)
    if mask.sum() > 5:
        zmid_run.append(10**((lo+hi)/2))
        frac_run.append(frac[mask].mean() * 100)
ax.plot(zmid_run, frac_run, color='darkorange', lw=2, label='Running mean')
ax.axhline(frac.mean() * 100, color='k', lw=1, ls='--',
           label=f'Global mean = {frac.mean()*100:.2f}%')
ax.set_xscale('log')
ax.set_xlabel(r'$z_\mathrm{CMB}$', fontsize=11)
ax.set_ylabel(r'$|\Delta\xi / \xi_\mathrm{mid}|$ (%)', fontsize=11)
ax.set_title(r'Fractional $\xi$-residual vs redshift (1,671 SNe Ia)', fontsize=11)
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(fig_path('xi_vs_z_scatter.png'), **STYLE)
plt.close()

# ----------------------------------------------------------------
# 3. CDF of |Δξ/ξ|
# ----------------------------------------------------------------
print("Fig 3: xi_residual_cdf ...")
fig, ax = plt.subplots(figsize=(6, 4))
sorted_frac = np.sort(frac * 100)
cdf = np.arange(1, len(sorted_frac)+1) / len(sorted_frac)
ax.plot(sorted_frac, cdf, color='steelblue', lw=2)
ax.axvline(frac.mean() * 100, color='darkorange', lw=1.5, ls='--',
           label=f'Mean = {frac.mean()*100:.2f}%')
ax.axhline(0.5, color='grey', lw=0.8, ls=':')
ax.axhline(0.9, color='grey', lw=0.8, ls=':')
p50 = sorted_frac[np.searchsorted(cdf, 0.5)]
p90 = sorted_frac[np.searchsorted(cdf, 0.9)]
ax.axvline(p50, color='grey', lw=0.8, ls=':', label=f'50th pct = {p50:.2f}%')
ax.axvline(p90, color='grey', lw=0.8, ls=':', label=f'90th pct = {p90:.2f}%')
ax.set_xlabel(r'$|\Delta\xi / \xi_\mathrm{mid}|$ (%)', fontsize=11)
ax.set_ylabel('CDF', fontsize=11)
ax.set_title(r'CDF of fractional $\xi$-residual (1,671 SNe Ia)', fontsize=11)
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(fig_path('xi_residual_cdf.png'), **STYLE)
plt.close()

# ----------------------------------------------------------------
# 4. chi2 sanity check (single Om)
# ----------------------------------------------------------------
print("Fig 4: chi2_sanity_check ...")
H0_grid = np.linspace(73.0, 78.5, 100)
chi2_grid = np.array([chi2_hf(h, Om_shoes, z_hf, mu_hf, sig_hf) for h in H0_grid])
chi2_min = chi2_grid.min()
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(H0_grid, chi2_grid - chi2_min, color='steelblue', lw=2)
ax.axvline(H0_s, color='k', lw=1.2, ls='--', label=f'Min: $H_0 = {H0_s:.2f}$')
ax.axhline(1.0, color='grey', lw=0.8, ls=':', label=r'$\Delta\chi^2 = 1$')
ax.set_xlabel(r'$H_0$ (km s$^{-1}$ Mpc$^{-1}$)', fontsize=11)
ax.set_ylabel(r'$\Delta\chi^2(H_0)$', fontsize=11)
ax.set_title(r'$\chi^2$ profile at $\Omega_m = 0.334$', fontsize=11)
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(fig_path('chi2_sanity_check.png'), **STYLE)
plt.close()

# ----------------------------------------------------------------
# 5. chi2 comparison (both Om values)
# ----------------------------------------------------------------
print("Fig 5: chi2_comparison ...")
chi2_p_grid = np.array([chi2_hf(h, Om_planck, z_hf, mu_hf, sig_hf) for h in H0_grid])
chi2_p_min = chi2_p_grid.min()
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(H0_grid, chi2_grid  - chi2_min,   color='steelblue',   lw=2,
        label=r'$\Omega_m = 0.334$')
ax.plot(H0_grid, chi2_p_grid - chi2_p_min, color='darkorange', lw=2, ls='--',
        label=r'$\Omega_m = 0.315$')
ax.axvline(H0_s, color='steelblue',  lw=1, ls=':')
ax.axvline(H0_p, color='darkorange', lw=1, ls=':')
ax.annotate('', xy=(H0_s, 0.5), xytext=(H0_p, 0.5),
            arrowprops=dict(arrowstyle='<->', color='k', lw=1.2))
ax.text((H0_s+H0_p)/2, 0.55,
        fr'$\Delta H_0 = {H0_s-H0_p:+.3f}$', ha='center', fontsize=9)
ax.set_xlabel(r'$H_0$ (km s$^{-1}$ Mpc$^{-1}$)', fontsize=11)
ax.set_ylabel(r'$\Delta\chi^2(H_0)$', fontsize=11)
ax.set_title(r'$\chi^2$ profiles at both $\Omega_m$ values', fontsize=11)
ax.legend(fontsize=9)
ax.set_ylim(-0.1, 3.0)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(fig_path('chi2_comparison.png'), **STYLE)
plt.close()

# ----------------------------------------------------------------
# 6. Hubble diagram residuals vs z
# ----------------------------------------------------------------
print("Fig 6: residual_vs_z ...")
mu_th = mu_theory(z_hf, H0_s, Om_shoes) + M_fixed
resid = mu_hf - mu_th
log_edges = np.linspace(np.log10(z_hf.min()), np.log10(z_hf.max()), 25)
zmid_run, res_run = [], []
for lo, hi in zip(log_edges[:-1], log_edges[1:]):
    mask = (np.log10(z_hf) >= lo) & (np.log10(z_hf) < hi)
    if mask.sum() > 5:
        zmid_run.append(10**((lo+hi)/2))
        res_run.append(resid[mask].mean())
fig, ax = plt.subplots(figsize=(7, 4))
ax.scatter(z_hf, resid, s=1, c='grey', alpha=0.3, rasterized=True)
ax.plot(zmid_run, res_run, color='darkorange', lw=2, label='Running mean')
ax.axhline(0, color='k', lw=0.8, ls='--')
ax.set_xscale('log')
ax.set_xlabel(r'$z_\mathrm{CMB}$', fontsize=11)
ax.set_ylabel(r'$\Delta\mu = \mu^\mathrm{obs} - \mu^\mathrm{th}$ (mag)', fontsize=11)
ax.set_title(r'Hubble diagram residuals at best-fit $H_0(\Omega_m=0.334)$', fontsize=11)
ax.legend(fontsize=9)
ax.set_ylim(-2.5, 2.5)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(fig_path('residual_vs_z.png'), **STYLE)
plt.close()

# ----------------------------------------------------------------
# 7. Calibrator CEPH_DIST vs MU_SH0ES
# ----------------------------------------------------------------
print("Fig 7: calibrator_ceph_vs_mu ...")
fig, ax = plt.subplots(figsize=(6, 5))
ax.errorbar(mu_cal, cal['MU_SH0ES'].values,
            xerr=sig_cal, fmt='o', color='steelblue',
            markersize=4, alpha=0.8, elinewidth=0.8, capsize=2)
lims = [min(mu_cal.min(), cal['MU_SH0ES'].values.min()) - 0.3,
        max(mu_cal.max(), cal['MU_SH0ES'].values.max()) + 0.3]
ax.plot(lims, lims, 'k--', lw=1, label='1:1')
ax.plot(lims, [x + M_fixed for x in lims], 'r-', lw=1.5,
        label=fr'$M_\mathrm{{fixed}} = {M_fixed:+.4f}$ mag')
ax.set_xlim(lims); ax.set_ylim(lims)
ax.set_xlabel(r'\texttt{CEPH\_DIST} (mag)', fontsize=11)
ax.set_ylabel(r'\texttt{MU\_SH0ES} (mag)', fontsize=11)
ax.set_title('Calibrator SNe Ia: Cepheid vs SH0ES distance modulus', fontsize=11)
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(fig_path('calibrator_ceph_vs_mu.png'), **STYLE)
plt.close()

# ----------------------------------------------------------------
# 8. z distribution of HF sample
# ----------------------------------------------------------------
print("Fig 8: z_distribution ...")
z_edges = np.percentile(z_hf, [0, 25, 50, 75, 100])
bin_colors = ['#4C72B0', '#DD8452', '#55A868', '#C44E52']
fig, ax = plt.subplots(figsize=(7, 4))
log_bins = np.logspace(np.log10(z_hf.min()*0.9), np.log10(z_hf.max()*1.1), 60)
for i, (lo, hi, c) in enumerate(zip(z_edges[:-1], z_edges[1:], bin_colors)):
    mask = (z_hf >= lo) & (z_hf < hi) if i < 3 else (z_hf >= lo)
    ax.hist(z_hf[mask], bins=log_bins, color=c, alpha=0.75,
            label=f'Bin {i+1}: z=[{lo:.3f},{hi:.3f}]')
ax.set_xscale('log')
ax.set_xlabel(r'$z_\mathrm{CMB}$', fontsize=11)
ax.set_ylabel('Count', fontsize=11)
ax.set_title('Redshift distribution of 1,616 Hubble-flow SNe Ia', fontsize=11)
ax.legend(fontsize=8, ncol=2)
ax.grid(alpha=0.3, axis='y')
fig.tight_layout()
fig.savefig(fig_path('z_distribution.png'), **STYLE)
plt.close()

# ----------------------------------------------------------------
# 9. σ_μ vs z
# ----------------------------------------------------------------
print("Fig 9: sigma_vs_z ...")
log_edges = np.linspace(np.log10(z_hf.min()), np.log10(z_hf.max()), 25)
zmid_run, sig_run = [], []
for lo, hi in zip(log_edges[:-1], log_edges[1:]):
    mask = (np.log10(z_hf) >= lo) & (np.log10(z_hf) < hi)
    if mask.sum() > 5:
        zmid_run.append(10**((lo+hi)/2))
        sig_run.append(sig_hf[mask].mean())
fig, ax = plt.subplots(figsize=(7, 4))
ax.scatter(z_hf, sig_hf, s=1, c='grey', alpha=0.3, rasterized=True)
ax.plot(zmid_run, sig_run, color='darkorange', lw=2, label='Running mean')
ax.set_xscale('log')
ax.set_xlabel(r'$z_\mathrm{CMB}$', fontsize=11)
ax.set_ylabel(r'$\sigma_\mu$ (mag)', fontsize=11)
ax.set_title(r'Distance-modulus uncertainty vs redshift (1,616 HF SNe Ia)', fontsize=11)
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(fig_path('sigma_vs_z.png'), **STYLE)
plt.close()

# ----------------------------------------------------------------
# 10. H0(Om) scan per z-bin
# ----------------------------------------------------------------
print("Fig 10: H0_scan_per_zbin ...")
Om_scan = np.linspace(0.28, 0.36, 17)
z_edges = np.percentile(z_hf, [0, 25, 50, 75, 100])
bin_colors = ['#4C72B0', '#DD8452', '#55A868', '#C44E52']
fig, ax = plt.subplots(figsize=(7, 5))
for i, (lo, hi, c) in enumerate(zip(z_edges[:-1], z_edges[1:], bin_colors)):
    mask = (z_hf >= lo) & (z_hf < hi) if i < 3 else (z_hf >= lo)
    zb, mub, sb = z_hf[mask], mu_hf[mask], sig_hf[mask]
    h0s = [best_H0(om, zb, mub, sb) for om in Om_scan]
    ax.plot(Om_scan, h0s, color=c, lw=1.8,
            label=f'Bin {i+1}: z=[{lo:.3f},{hi:.3f}]')
ax.axvline(Om_shoes,  color='blue',   lw=1, ls='--', alpha=0.7, label=r'SH0ES $\Omega_m=0.334$')
ax.axvline(Om_planck, color='orange', lw=1, ls='--', alpha=0.7, label=r'Planck $\Omega_m=0.315$')
ax.set_xlabel(r'$\Omega_m$', fontsize=11)
ax.set_ylabel(r'$H_0$ (km s$^{-1}$ Mpc$^{-1}$)', fontsize=11)
ax.set_title(r'$H_0(\Omega_m)$ scan per redshift quartile', fontsize=11)
ax.legend(fontsize=8, ncol=2)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(fig_path('H0_scan_per_zbin.png'), **STYLE)
plt.close()

# ----------------------------------------------------------------
# 11. H0 comparison ladder
# ----------------------------------------------------------------
print("Fig 11: H0_comparison_ladder ...")
# Literature values: (label, H0, err_lo, err_hi, colour)
literature = [
    ('Planck 2018\n(CMB)',              67.4, 0.5,  0.5,  'royalblue'),
    ('WMAP 9-yr\n(CMB)',                69.32, 0.80, 0.80, 'cornflowerblue'),
    ('Freedman+2021\n(TRGB)',           69.8, 1.7,  1.7,  'mediumseagreen'),
    ('Freedman+2024\n(TRGB/JWST)',      70.7, 2.1,  2.1,  'mediumseagreen'),
    ('Birrer+2020\n(TDCOSMO)',          74.5, 5.6,  5.6,  'mediumorchid'),
    ('Wong+2020\n(H0LiCOW)',            73.3, 1.7,  1.8,  'mediumorchid'),
    ('Pesce+2020\n(Megamasers)',        73.9, 3.0,  3.0,  'peru'),
    ('Blakeslee+2021\n(SBF)',           73.3, 2.1,  2.1,  'sandybrown'),
    ('Riess+2022\n(SH0ES)',             73.04, 1.04, 1.04, 'crimson'),
]
fig, ax = plt.subplots(figsize=(7, 6))
y = np.arange(len(literature), 0, -1)
for yi, (label, h0, elo, ehi, c) in zip(y, literature):
    ax.errorbar(h0, yi, xerr=[[elo], [ehi]], fmt='o', color=c,
                markersize=6, elinewidth=2, capsize=4)
ax.yaxis.set_ticks(y)
ax.yaxis.set_ticklabels([l for l, *_ in literature], fontsize=8)
ax.axvspan(H0_s - 0.05, H0_s + 0.05, color='steelblue', alpha=0.25,
           label=fr'This work: $H_0 = {H0_s:.2f}$ ($\Omega_m=0.334$)')
ax.axvspan(H0_p - 0.05, H0_p + 0.05, color='darkorange', alpha=0.25,
           label=fr'This work: $H_0 = {H0_p:.2f}$ ($\Omega_m=0.315$)')
ax.set_xlabel(r'$H_0$ (km s$^{-1}$ Mpc$^{-1}$)', fontsize=11)
ax.set_title(r'Published $H_0$ measurements (1$\sigma$ error bars)', fontsize=11)
ax.legend(fontsize=8, loc='lower right')
ax.grid(alpha=0.3, axis='x')
fig.tight_layout()
fig.savefig(fig_path('H0_comparison_ladder.png'), **STYLE)
plt.close()

# ----------------------------------------------------------------
# 12. xi residual stats per z-bin
# ----------------------------------------------------------------
print("Fig 12: ximid_vs_z_bins ...")
z_edges_full = np.percentile(z_full, [0, 25, 50, 75, 100])
bin_fracs = []
for i, (lo, hi) in enumerate(zip(z_edges_full[:-1], z_edges_full[1:])):
    mask = (z_full >= lo) & (z_full < hi) if i < 3 else (z_full >= lo)
    bin_fracs.append(frac[mask].mean() * 100)
fig, ax = plt.subplots(figsize=(7, 4))
colors = ['#4C72B0', '#DD8452', '#55A868', '#C44E52']
labels = [f'z=[{z_edges_full[i]:.3f},{z_edges_full[i+1]:.3f}]' for i in range(4)]
ax.bar(range(1, 5), bin_fracs, color=colors, alpha=0.85, edgecolor='white')
ax.axhline(frac.mean() * 100, color='k', lw=1.5, ls='--',
           label=f'Global mean = {frac.mean()*100:.2f}%')
ax.set_xticks(range(1, 5))
ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel(r'$\langle|\Delta\xi|\rangle / \langle\xi\rangle$ (%)', fontsize=11)
ax.set_title(r'$\xi$-space residual per redshift quartile', fontsize=11)
ax.legend(fontsize=9)
ax.grid(alpha=0.3, axis='y')
fig.tight_layout()
fig.savefig(fig_path('ximid_vs_z_bins.png'), **STYLE)
plt.close()

print("\nAll 12 figures generated.")
print(f"Saved to: {HERE}")
