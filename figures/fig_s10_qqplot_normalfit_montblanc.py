"""Plotting of Figure S10: Q-Q plot and normal fit after standardization of dh for the Mont-Blanc case study"""
import numpy as np
import xdem
import geoutils as gu
from scipy.stats import probplot, norm
import matplotlib.pyplot as plt

fn_ddem = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/dh_Pleiades-SPOT6_Mont-Blanc_NK_Deramp.tif'
fn_pleiades = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/Pleiades_Mont-Blanc_2017-10-25_DEM_5m.tif'
fn_shp = '/home/atom/data/inventory_products/RGI/00_rgi60_neighb_merged/11_rgi60_CentralEurope/11_rgi60_CentralEurope.shp'
fn_forest = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/outlines/forest_Mont-Blanc_ESACCI_delainey.shp'

pleia_ddem = gu.Raster(fn_ddem)
ref_dem = gu.Raster(fn_pleiades)
glaciers_outlines = gu.Vector(fn_shp)
forest_outlines = gu.Vector(fn_forest)
mask_glacier = glaciers_outlines.create_mask(pleia_ddem)
mask_forest = forest_outlines.create_mask(pleia_ddem)

# Remove forest, very large outliers
pleia_ddem.data[mask_forest] = np.nan
pleia_ddem.data[np.abs(pleia_ddem.data)>200] = np.nan

# pleia_ddem.data[mask_glacier] = np.nan

slope, planc, profc = xdem.terrain.get_terrain_attribute(ref_dem, attribute=['slope', 'planform_curvature',
                                                                                'profile_curvature'])
maxabsc = np.maximum(np.abs(planc), np.abs(profc))

del ref_dem

# # Filter large outliers per category
bins_slope = [0, 2.5, 5, 10, 15, 20, 30, 40, 50, 70, 90]
# bins_curv = [-10, -8, -6, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 6, 8, 10]
bins_curv = [0, 0.1, 0.2, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 20, 50]
for i in range(len(bins_slope) - 1):
    # Subset by slope category
    subset = np.logical_and(slope.data >= bins_slope[i], slope.data < bins_slope[i + 1])
    dh_sub = pleia_ddem.data[subset]
    # Remove very large outliers of the category
    med_sub = np.nanmedian(dh_sub)
    nmad_sub = xdem.spatialstats.nmad(dh_sub)
    # Remove outliers
    pleia_ddem.data[np.logical_and(subset, np.abs(pleia_ddem.data-med_sub) > 5 * nmad_sub)] = np.nan
for i in range(len(bins_curv) - 1):
    # Subset by slope category
    subset = np.logical_and(maxabsc >= bins_curv[i], maxabsc < bins_curv[i + 1])
    dh_sub = pleia_ddem.data[subset]
    # Remove very large outliers of the category
    med_sub = np.nanmedian(dh_sub)
    nmad_sub = xdem.spatialstats.nmad(dh_sub)
    # Remove outliers
    pleia_ddem.data[np.logical_and(subset, np.abs(pleia_ddem.data-med_sub) > 5 * nmad_sub)] = np.nan


pleia_ddem_sta = pleia_ddem.data[~mask_glacier]
slope_sta = slope.data[~mask_glacier]
maxabsc_sta = maxabsc[~mask_glacier]

df_sub = xdem.spatialstats.nd_binning(pleia_ddem_sta, list_var=[slope_sta, maxabsc_sta], list_var_names=['slope', 'maxc'], list_var_bins=(bins_slope, bins_curv))

fn = xdem.spatialstats.interp_nd_binning(df_sub, list_var_names=['slope', 'maxc'], statistic='nmad', min_count=30)

maxabsc[maxabsc>50] = 50
dh_err = fn((slope.data, maxabsc))

pleia_ddem_gla = pleia_ddem.data[mask_glacier]

std_dh = pleia_ddem.data.data/dh_err
std_dh[np.abs(std_dh- np.nanmedian(std_dh))>5*xdem.spatialstats.nmad(std_dh)] = np.nan

std_dh_sta = std_dh[~mask_glacier]
std_dh_gla = std_dh[mask_glacier]

subsample_dh_sta = gu.spatial_tools.subsample_raster(pleia_ddem_sta, subsample=500000, random_state=42)
subsample_dh_gla = gu.spatial_tools.subsample_raster(pleia_ddem_gla, subsample=500000, random_state=42)

subsample_std_sta = gu.spatial_tools.subsample_raster(std_dh_sta, subsample=500000, random_state=42)
subsample_std_gla = gu.spatial_tools.subsample_raster(std_dh_gla, subsample=500000, random_state=42)


# 1/ First, we plot the Q-Q plot for elevation differences
fig = plt.figure(figsize=(7, 7))

grid = plt.GridSpec(22, 23, wspace=0.1, hspace=0.1)

ax = fig.add_subplot(grid[0:10, 0:10])

probplot(x=subsample_dh_sta, plot=plt, fit=False)

probplot(x=subsample_dh_gla, plot=plt, fit=False)

ax.get_lines()[0].set_markerfacecolor('tab:brown')
ax.get_lines()[0].set_markeredgecolor('tab:brown')
ax.get_lines()[0].set_markersize(4.0)

# ax.get_lines()[1].set_color('black')

ax.get_lines()[1].set_markerfacecolor('tab:cyan')
ax.get_lines()[1].set_markeredgecolor('tab:cyan')
ax.get_lines()[1].set_markersize(4.0)

ax.plot([-5, 5], [-5, 5], color='black', zorder=1)

ax.set_title('')
ax.set_xlabel('Theoretical normal quantiles')
ax.set_ylabel('Quantiles of ordered\nelevation differences')

ax.plot([], [], color='tab:brown', label='Stable terrain', linewidth=5)
ax.plot([], [], color='tab:cyan', label='Moving terrain', linewidth=5)
ax.plot([], [], color='black', label='1:1 line')

ax.legend(loc='lower right')

ax.text(0.05, 0.95, 'a', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14,
         zorder=20)

# 2/ Then with the standardized elevation differences
ax = fig.add_subplot(grid[12:, 0:10])

probplot(x=subsample_std_sta, plot=plt, fit=False)

probplot(x=subsample_std_gla, plot=plt, fit=False)

ax.get_lines()[0].set_markerfacecolor('tab:brown')
ax.get_lines()[0].set_markeredgecolor('tab:brown')
ax.get_lines()[0].set_markersize(4.0)

# ax.get_lines()[1].set_color('black')

ax.get_lines()[1].set_markerfacecolor('tab:cyan')
ax.get_lines()[1].set_markeredgecolor('tab:cyan')
ax.get_lines()[1].set_markersize(4.0)
ax.set_xlabel('Theoretical normal quantiles')
ax.set_ylabel('Quantiles of ordered\nstandard score of\nelevation differences')

ax.plot([-5, 5], [-5, 5], color='black', zorder=1)

ax.set_title('')
ax.text(0.05, 0.95, 'c', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14,
         zorder=20)

# 3/ Histogram comparison with a normal fit for the elevation differences
ax = fig.add_subplot(grid[0:10, 13:18])

n, bins, patches = ax.hist(subsample_dh_sta, color='tab:brown', alpha=0.6, bins=50, density=True, range=(-20, 20))
ax.hist(subsample_dh_gla, color='tab:cyan', alpha=0.6, bins=50, density=True, range=(-20, 20))

mu_dh, sigma_dh = norm.fit(subsample_dh_sta)

y = norm.pdf(bins, mu_dh, sigma_dh)
ax.plot(bins, y, 'black', linewidth=1, linestyle='dashed', label='Normal\nfit')
ax.text(0.1, 0.95, 'b', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14,
         zorder=20)
ax.set_xlim((-20, 0))
ax.set_ylabel('Probability density')
ax.set_xlabel('Elevation differences (m)', x=1, ha='center')
ax.legend(loc='center left')

ax = fig.add_subplot(grid[0:10, 18:])

ax.hist(subsample_dh_sta, color='tab:brown', alpha=0.6, bins=50, density=True, range=(-20, 20))
ax.hist(subsample_dh_gla, color='tab:cyan', alpha=0.6, bins=50, density=True, range=(-20, 20))
ax.plot(bins, y, 'black', linewidth=1, linestyle='dashed')

ax.set_xlim((0.01, 20))
ax.set_yscale('log')
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")

# 4/ Same for the standardized differences
ax = fig.add_subplot(grid[12:, 13:18])

n, bins, patches = ax.hist(subsample_std_sta, color='tab:brown', alpha=0.6, bins=50, density=True)
ax.hist(subsample_std_gla, color='tab:cyan', alpha=0.6, bins=50, density=True)

mu_std, sigma_std = norm.fit(subsample_std_sta)

y = norm.pdf(bins, mu_std, sigma_std)
ax.plot(bins, y, 'black', linewidth=1, linestyle='dashed')
ax.text(0.1, 0.95, 'd', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14,
         zorder=20)
# ax.set_yscale('log')
ax.set_xlim((-5, 0))
ax.set_ylabel('Probability density')
ax.set_xlabel('Standardized elevation differences', x=1, ha='center')

ax = fig.add_subplot(grid[12:, 18:])

ax.hist(subsample_std_sta, color='tab:brown', alpha=0.6, bins=50, density=True, range=(-20, 20))
ax.hist(subsample_std_gla, color='tab:cyan', alpha=0.6, bins=50, density=True, range=(-20, 20))
ax.plot(bins, y, 'black', linewidth=1, linestyle='dashed')

ax.set_xlim((0.01, 5))
ax.set_yscale('log')
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")

# Save to file
plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_S10_final.png', dpi=400)

