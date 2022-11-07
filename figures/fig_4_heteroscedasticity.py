"""Plotting of Figure 4: elevation heteroscedasticity for the Mont-Blanc case study"""
import numpy as np
import pandas as pd
import xdem
import geoutils as gu
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import seaborn as sns
from matplotlib.collections import LineCollection
from matplotlib.legend_handler import HandlerLineCollection, HandlerPatch
import matplotlib.colors as colors
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature

# No need to run any estimation script for this plotting function, that redoes the same processing
# (this is because the violin plots need all samples, not just a binned estimation of heteroscedasticity)

# Open files
fn_ddem = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/dh_Pleiades-SPOT6_Mont-Blanc_NK_Deramp.tif'
fn_pleiades = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/Pleiades_Mont-Blanc_2017-10-25_DEM_5m.tif'
fn_shp = '/home/atom/data/inventory_products/RGI/00_rgi60_neighb_merged/11_rgi60_CentralEurope/11_rgi60_CentralEurope.shp'
fn_forest = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/outlines/forest_Mont-Blanc_ESACCI_delainey.shp'
fn_hs = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/Mont-Blanc_2017-10-25_DEM_5m_hillshade.tif'

pleia_ddem = gu.Raster(fn_ddem)
ref_dem = gu.Raster(fn_pleiades)
glaciers_outlines = gu.Vector(fn_shp)
forest_outlines = gu.Vector(fn_forest)
mask_glacier = glaciers_outlines.create_mask(pleia_ddem)
mask_forest = forest_outlines.create_mask(pleia_ddem)

# Remove forest, very large outliers
pleia_ddem.data[mask_forest] = np.nan
pleia_ddem.data[np.abs(pleia_ddem.data)>500] = np.nan

# Derive terrain slope, planform and profile curvature, then maximum curvature
slope, planc, profc = xdem.terrain.get_terrain_attribute(ref_dem, attribute=['slope', 'planform_curvature',
                                                                                'profile_curvature'])
maxabsc = np.maximum(np.abs(planc), np.abs(profc))

# Define bins, and filter very large outliers per category
bins_slope = [0, 2.5, 5, 10, 15, 20, 30, 40, 50, 70, 90]
bins_curv = [0, 0.2, 0.5, 1, 2, 3, 4, 6, 10, 20, 50]
for i in range(len(bins_slope) - 1):
    # Subset by slope category
    subset = np.logical_and(slope.data >= bins_slope[i], slope.data < bins_slope[i + 1])
    dh_sub = pleia_ddem.data[subset]
    # Remove very large outliers of the category
    med_sub = np.nanmedian(dh_sub)
    nmad_sub = xdem.spatialstats.nmad(dh_sub)
    pleia_ddem.data[np.logical_and(subset, np.abs(pleia_ddem.data-med_sub) > 5 * nmad_sub)] = np.nan
for i in range(len(bins_curv) - 1):
    # Subset by curvature category
    subset = np.logical_and(maxabsc >= bins_curv[i], maxabsc < bins_curv[i + 1])
    dh_sub = pleia_ddem.data[subset]
    # Remove very large outliers of the category
    med_sub = np.nanmedian(dh_sub)
    nmad_sub = xdem.spatialstats.nmad(dh_sub)
    pleia_ddem.data[np.logical_and(subset, np.abs(pleia_ddem.data-med_sub) > 5 * nmad_sub)] = np.nan

# Bin in 1D with slope
df = xdem.spatialstats.nd_binning(pleia_ddem.data.ravel(), list_var=[slope.data.ravel()], list_var_names=['slope'])

# Format slope data for seaborn plots (will require all samples becaue of the violin plot)
terrain = np.array(['Moving' if mask_glacier.ravel()[i] else 'Stable' for i in range(len(mask_glacier.ravel()))])
list_df, list_counts, list_mean, list_std, list_bin_name_slope, list_bin_name_2_slope = ([] for i in range(6))
for i in range(len(bins_slope)-1):
    # Subset by slope category
    subset = np.logical_and(slope.data.ravel() >= bins_slope[i], slope.data.ravel() < bins_slope[i+1])
    dh_sub = pleia_ddem.data.ravel()[subset]
    terrain_sub = terrain[subset.ravel()]
    # Remove very large outliers of the category
    dh_sub[np.abs(dh_sub - np.nanmedian(dh_sub))>5*xdem.spatialstats.nmad(dh_sub)] = np.nan
    # Store in dataframe
    df_subset = pd.DataFrame()
    df_subset = df_subset.assign(dh = dh_sub, terrain = terrain_sub)
    bin_name = str(bins_slope[i])+'–\n'+str(bins_slope[i+1])
    bin_name_2 = str(bins_slope[i])+'–'+str(bins_slope[i+1])
    df_subset['slope_category'] = bin_name
    list_df.append(df_subset)
    # Store counts for histograms
    count_unstable = np.sum(np.isfinite(dh_sub[terrain_sub=='Moving']))
    count_stable = np.sum(np.isfinite(dh_sub[terrain_sub=='Stable']))
    list_counts.append((count_unstable, count_stable))
    mean = np.nanmean(dh_sub)
    std = xdem.spatialstats.nmad(dh_sub)
    list_mean.append(mean)
    list_std.append(std)
    list_bin_name_slope.append(bin_name)
    list_bin_name_2_slope.append(bin_name_2)

list_mean = np.array(list_mean)
list_std = np.array(list_std)

# Concatenate all dataframes
df_all = pd.concat(list_df)

# Prepare a consistent colormap for the figure
cmap = plt.get_cmap('YlOrRd', 100)
# cmap._init()
# cmap._lut[:, -1] = 1.
col_bounds = np.array([0.8, 1.2, 2, 3.5, 5])
cb = []
cb_val = np.linspace(0, 1, len(col_bounds))
for j in range(len(cb_val)):
    cb.append(cmap(cb_val[j]))
cmap_cus = colors.LinearSegmentedColormap.from_list('my_cb', list(
    zip((col_bounds - min(col_bounds)) / (max(col_bounds - min(col_bounds))), cb)), N=1000)

# 0/ Start figure
fig = plt.figure(figsize=(15, 5.25))
grid = plt.GridSpec(24, 32, wspace=0.1, hspace=0.1)

# 1/ Panel a on slope

# 1.1/ First, an horizontal axis on top to plot the sample histograms
ax0 = fig.add_subplot(grid[:2, :10])

for i in range(len(list_counts)):
    c = list_counts[i]
    c_unstable = c[0]
    c_stable = c[1]
    ax0.fill_between([i, i+0.3], [0] * 2, [c_unstable] * 2,
                     facecolor='tab:cyan', edgecolor='black', alpha=1,
                     linewidth=0.5)
    ax0.fill_between([i-0.3, i], [0] * 2, [c_stable] * 2,
                     facecolor='tab:brown', edgecolor='black', alpha=1,
                     linewidth=0.5)

ax0.set_ylabel('Sample\ncount')
ax0.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.set_xlim((-0.5, len(list_counts)-0.5))
ax0.set_yscale('log')
# ax0.set_ylim((10, np.nanmax(np.array(list_counts))))
ax0.set_xticks([])
ax0.text(0.035, 0.99, 'a', transform=ax0.transAxes, ha='left', va='top',fontweight='bold',fontsize=14)

# 1.2/ Then, another axis to plot the violin plots
ax = fig.add_subplot(grid[2:9, :10])


# Plot a colored curve for the dispersion, get a linearized version of the bins to convert from discrete to
# continuous plotting
x0 = np.arange(0, len(list_counts))
bins_slope_center = np.subtract(np.asarray(bins_slope[1:]), np.diff(bins_slope) / 2)
interp_x_y = interp1d(x0, bins_slope_center)
xnew = np.linspace(0, len(list_counts)-1, 1000)
ynew = interp_x_y(xnew)

# Fit an exponential function to the dispersion
def fun(x, a, b):
    return a + np.sqrt(b*np.tan(np.deg2rad(x)))

res, _ = curve_fit(fun, xdata=bins_slope_center, ydata=list_std)

ax.hlines(0, xmin=-1, xmax=15, color='tab:grey', linestyles='dashed', zorder=1)

# yfit = fun(ynew, *res)
piecewise = interp1d(x0, list_std, fill_value='extrapolate')
yfit = piecewise(xnew)
points = np.array([xnew, 2*yfit]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
norm = plt.Normalize(0, 10)
lc = LineCollection(segments, cmap=cmap_cus, zorder=2)
lc.set_array(yfit)
lc.set_linewidth(5)
ax.add_collection(lc)

a = list(ax.get_children())

# Finally, plot the violin plot
sns.violinplot(ax=ax, x="slope_category", y="dh", hue="terrain",
               data=df_all, split=True,  palette={'Moving':'tab:cyan', 'Stable':'tab:brown'}, linewidth=1,
               edgecolor='black')

ax.set_ylim((-15, 15))

# Clear slope dataset to limit memory usage
del df_all

b = list(ax.get_children())
c = set(b) - set(a)
for d in c:
    d.set_zorder(1000)

ax.set_xlabel('Slope categories (degrees)')
ax.set_ylabel('Elevation\ndifferences\n(m)', labelpad=-5)
ax.set_xlim((-0.5, len(list_counts)-0.5))

# Customize legend
p0 = Patch(facecolor='tab:brown', edgecolor='black', label='Stable terrain')
p1 = Patch(facecolor='tab:cyan', edgecolor='black', label='Moving terrain')

class HandlerColorLineCollection(HandlerLineCollection):
    def create_artists(self, legend, artist, xdescent, ydescent, width, height, fontsize, trans):
        x = np.linspace(0, width, self.get_numpoints(legend) + 1)
        y = np.zeros(self.get_numpoints(legend) + 1) + height / 2. - ydescent
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[: -1], points[1: ]], axis = 1)
        lc = LineCollection(segments, cmap = artist.cmap,
        transform = trans)
        lc.set_array(x)
        lc.set_linewidth(artist.get_linewidth())
        return [lc]

points = np.array([[0, 0.5, 1], [0, 0, 0]]).T.reshape(-1, 1, 2)
segments=np.concatenate([points[:-1], points[1:]], axis=1)
lc_leg = LineCollection(segments=segments, cmap=cmap_cus)
lc_leg.set_array(np.array([0, 5, 10]))
lc_leg.set_linewidth(5)

l0 = ax.legend([lc_leg], ['Dispersion (2$\sigma$)'],
           handler_map={lc_leg: HandlerColorLineCollection(numpoints=100)},
           framealpha=1, loc='upper left', ncol=1)
l0.set_zorder(30)

l1 = ax.legend([p0, p1], ['Stable', 'Moving'],
           handler_map={p0 : HandlerPatch(), p1: HandlerPatch()},
           framealpha=1, loc='lower left', ncol=2)
l1.set_zorder(30)
ax.add_artist(l0)

# Format curvature data for seaborn plots (will require all samples becaue of the violin plot)
terrain = np.array(['Moving' if mask_glacier.ravel()[i] else 'Stable' for i in range(len(mask_glacier.ravel()))])
list_df, list_counts, list_mean, list_std, list_bin_name_curv = ([] for i in range(5))
for i in range(len(bins_curv)-1):
    # Subset by curvature category
    subset = np.logical_and(maxabsc.ravel() >= bins_curv[i], maxabsc.ravel() < bins_curv[i+1])
    dh_sub = pleia_ddem.data.ravel()[subset]
    terrain_sub = terrain[subset.ravel()]
    # Remove very large outliers of the category
    dh_sub[np.abs(dh_sub - np.nanmedian(dh_sub))>5*xdem.spatialstats.nmad(dh_sub)] = np.nan
    # Store in dataframe
    df_subset = pd.DataFrame()
    df_subset = df_subset.assign(dh = dh_sub, terrain = terrain_sub)
    if i <= 2:
        bin_name =  '{:.1f}'.format(bins_curv[i])+'–\n'+'{:.1f}'.format(bins_curv[i+1])
    elif len(bins_curv)-i-1<=3:
        bin_name =  '{:.0f}'.format(bins_curv[i])+'–\n'+'{:.0f}'.format(bins_curv[i+1])
    else:
        bin_name =  '{:.0f}'.format(bins_curv[i])+'–'+'{:.0f}'.format(bins_curv[i+1])


    df_subset['curv_category'] = bin_name
    list_df.append(df_subset)
    # Store counts for histograms
    count_unstable = np.sum(np.isfinite(dh_sub[terrain_sub=='Moving']))
    count_stable = np.sum(np.isfinite(dh_sub[terrain_sub=='Stable']))
    list_counts.append((count_unstable, count_stable))
    mean = np.nanmean(dh_sub)
    std = xdem.spatialstats.nmad(dh_sub)
    list_mean.append(mean)
    list_std.append(std)
    list_bin_name_curv.append(bin_name)

list_mean = np.array(list_mean)
list_std = np.array(list_std)

# Concatenate all dataframes
df_all_curv = pd.concat(list_df)

# 2/ Panel b on curvature
# 2.2/ First, an horizontal axis on top to plot the sample histograms
ax0 = fig.add_subplot(grid[13:15, :10])

for i in range(len(list_counts)):
    c = list_counts[i]
    c_unstable = c[0]
    c_stable = c[1]
    ax0.fill_between([i, i+0.3], [0] * 2, [c_unstable] * 2,
                     facecolor='tab:cyan', edgecolor='black', alpha=1,
                     linewidth=0.5)
    ax0.fill_between([i-0.3, i], [0] * 2, [c_stable] * 2,
                     facecolor='tab:brown', edgecolor='black', alpha=1,
                     linewidth=0.5)

ax0.set_ylabel('Sample\ncount')
ax0.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.set_xlim((-0.5, len(list_counts)-0.5))
ax0.set_yscale('log')
# ax0.set_ylim((10, np.nanmax(np.array(list_counts))))
ax0.set_xticks([])
ax0.text(0.035, 0.99, 'b', transform=ax0.transAxes, ha='left', va='top',fontweight='bold',fontsize=14)

# 2.2/ Then, another axis to plot the violin plots
ax = fig.add_subplot(grid[15:22, :10])


# Get a linearized version of the bins to convert from discrete to continuous plotting
x0 = np.arange(0, len(list_counts))
bins_curv_center = np.subtract(np.asarray(bins_curv[1:]), np.diff(bins_curv) / 2)
interp_x_y = interp1d(x0, bins_curv_center)
xnew = np.linspace(0, len(list_counts)-1, 1000)
ynew = interp_x_y(xnew)

# Fit a linear function to the dispersion
def fun(x, a, b):
    return a*x + b

res, _ = curve_fit(fun, xdata=bins_curv_center, ydata=list_std)

ax.hlines(0, xmin=-1, xmax=15, color='tab:grey', linestyles='dashed', zorder=1)

piecewise = interp1d(x0, list_std, fill_value='extrapolate')
yfit = piecewise(xnew)
points = np.array([xnew, 2*yfit]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
lc = LineCollection(segments, cmap=cmap_cus, zorder=2)
lc.set_array(yfit)
lc.set_linewidth(5)
ax.add_collection(lc)

a = list(ax.get_children())

sns.violinplot(ax=ax, x="curv_category", y="dh", hue="terrain",
               data=df_all_curv, split=True,  palette={'Moving':'tab:cyan', 'Stable':'tab:brown'}, linewidth=1,
               edgecolor='black')

ax.set_ylim((-15, 15))

# Clear curvature dataset to limit memory usage
del df_all_curv

b = list(ax.get_children())
c = set(b)-set(a)
for d in c:
    d.set_zorder(1000)

ax.set_xlabel('Maximum curvature categories (10$^{2}$ m$^{-1}$)')
ax.set_ylabel('Elevation\ndifferences\n(m)',labelpad=-5)
leg = ax.legend(loc='lower left')
leg.remove()
ax.set_xlim((-0.5, len(list_counts)-0.5))

# 3/ Plot panel c: the carpet with double distribution with slope and curvature
pleia_ddem_sta = pleia_ddem.data[~mask_glacier]
slope_sta = slope.data[~mask_glacier]
maxabsc_sta = maxabsc[~mask_glacier]

pleia_ddem_gla = pleia_ddem.data[mask_glacier]
slope_gla = slope.data[mask_glacier]
maxabsc_gla = maxabsc[mask_glacier]

df_sub_all = xdem.spatialstats.nd_binning(pleia_ddem.data, list_var=[slope.data, maxabsc], list_var_names=['slope', 'maxc'], list_var_bins=(bins_slope, bins_curv))
df_sub_all.loc[df_sub_all['count'].values < 30, ['nmad']] = np.nan
df_sub = xdem.spatialstats.nd_binning(pleia_ddem_sta, list_var=[slope_sta, maxabsc_sta], list_var_names=['slope', 'maxc'], list_var_bins=(bins_slope, bins_curv))
df_sub.loc[df_sub['count'].values < 30, ['nmad']] = np.nan
df_sub['slope_mid'] = pd.IntervalIndex(df_sub.slope).mid.values
df_sub['maxc_mid'] = pd.IntervalIndex(df_sub.maxc).mid.values
# df_sub.to_csv('/home/atom/ongoing/work_stderr_dem/case_study_montblanc/df_heteroscedas_slope_curv.csv')

df_sub_gla = xdem.spatialstats.nd_binning(pleia_ddem_gla, list_var=[slope_gla, maxabsc_gla], list_var_names=['slope', 'maxc'], list_var_bins=(bins_slope, bins_curv))
df_sub_gla.loc[df_sub_gla['count'].values < 30, 'nmad'] = np.nan


import warnings

scale_var_1 = 'linear'
scale_var_2 = 'linear'

var_name_1 = 'slope'
var_name_2 = 'maxc'

vmin = 0
vmax = 10

cmap = plt.get_cmap('YlOrRd')
statistic_name = 'nmad'

label_var_name_1 = 'Slope categories (degrees)'
label_var_name_2 = 'Maximum curvature categories (10$^{2}$ m$^{-1}$)'

label_statistic = 'Dispersion (1$\sigma$) of\nelevation\ndifferences (m)'

nodata_color = 'white'

# 3.1/ First, an horizontal axis on top to plot the sample histogram of the first variable
ax0 = fig.add_subplot(grid[:3, 12:19])
ax0.set_xscale(scale_var_1)
ax0.set_xticklabels([])
ax0.set_yscale('log')

# Plot the histogram manually with fill_between
interval_var_1 = pd.IntervalIndex(df_sub[var_name_1])
interval_var_1_gla = pd.IntervalIndex(df_sub_gla[var_name_1])
interval_var_1_all = pd.IntervalIndex(df_sub_all[var_name_1])
df_sub['var1_mid'] = interval_var_1.mid.values
df_sub_gla['var1_mid'] = interval_var_1_gla.mid.values
df_sub_all['var1_mid'] = interval_var_1_all.mid.values
unique_var_1 = np.unique(df_sub.var1_mid[np.isfinite(df_sub.var1_mid)])
list_counts = []
for i in range(len(unique_var_1)):
    df_var1 = df_sub[df_sub.var1_mid == unique_var_1[i]]
    df_var1_gla = df_sub_gla[df_sub_gla.var1_mid == unique_var_1[i]]
    count = np.nansum(df_var1['count'].values)
    count_gla = np.nansum(df_var1_gla['count'].values)
    list_counts.append(count)
    list_counts.append(count_gla)
    ax0.fill_between([i, i+0.3], [0] * 2, [count_gla] * 2,
                     facecolor='tab:cyan', edgecolor='black', alpha=1,
                     linewidth=0.5)
    ax0.fill_between([i-0.3, i], [0] * 2, [count] * 2,
                     facecolor='tab:brown', edgecolor='black', alpha=1,
                     linewidth=0.5)

ax0.set_ylabel('Sample\ncount')
# In case the axis value does not agree with the scale (e.g., 0 for log scale)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ax0.set_ylim((0, 1.1 * np.max(list_counts)))
    ax0.set_xlim((-0.5, len(bins_slope)-1.5))
# ax0.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
# Try to identify if the count is always the same
if np.sum(~(np.abs(list_counts[0] - np.array(list_counts)) < 5)) <= 2:
    ax0.text(0.5, 0.5, "Fixed number of\nsamples: " + '{:,}'.format(int(list_counts[0])), ha='center', va='center',
             fontweight='bold', transform=ax0.transAxes, bbox=dict(facecolor='white', alpha=0.8))

ax0.text(0.05, 0.90, 'c', transform=ax0.transAxes, ha='left', va='top',fontweight='bold',fontsize=14)


# 3.2/ Second, a vertical axis on the right to plot the sample histogram of the second variable
ax1 = fig.add_subplot(grid[3:22, 19:21])
ax1.set_yscale(scale_var_2)
ax1.set_yticklabels([])
ax1.set_xscale('log')

# Plot the histogram manually with fill_between
interval_var_2 = pd.IntervalIndex(df_sub[var_name_2])
interval_var_2_gla = pd.IntervalIndex(df_sub_gla[var_name_2])
interval_var_2_all = pd.IntervalIndex(df_sub_all[var_name_2])
df_sub['var2_mid'] = interval_var_2.mid.values
df_sub_gla['var2_mid'] = interval_var_2_gla.mid.values
df_sub_all['var2_mid'] = interval_var_2_all.mid.values
unique_var_2 = np.unique(df_sub.var2_mid[np.isfinite(df_sub.var2_mid)])
list_counts = []
for i in range(len(unique_var_2)):
    df_var2 = df_sub[df_sub.var2_mid == unique_var_2[i]]
    df_var2_gla = df_sub_gla[df_sub_gla.var2_mid == unique_var_2[i]]
    count = np.nansum(df_var2['count'].values)
    count_gla = np.nansum(df_var2_gla['count'].values)
    list_counts.append(count)
    list_counts.append(count_gla)
    ax1.fill_between([0, count_gla], [i] * 2, [i+0.3] * 2,
                     facecolor='tab:cyan', edgecolor='black', alpha=1,
                     linewidth=0.5)
    ax1.fill_between([0, count], [i-0.3] * 2, [i] * 2,
                     facecolor='tab:brown', edgecolor='black', alpha=1,
                     linewidth=0.5)

ax1.set_xlabel('Sample\ncount')
# In case the axis value does not agree with the scale (e.g., 0 for log scale)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ax1.set_xlim((0, 1.1 * np.max(list_counts)))
    ax1.set_ylim((-0.5, len(bins_curv)-1.5))
# ax1.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
# Try to identify if the count is always the same
if np.sum(~(np.abs(list_counts[0] - np.array(list_counts)) < 5)) <= 2:
    ax1.text(0.5, 0.5, "Fixed number of\nsamples: " + '{:,}'.format(int(list_counts[0])), ha='center', va='center',
             fontweight='bold', transform=ax1.transAxes, rotation=90, bbox=dict(facecolor='white', alpha=0.8))

# 3.3/ Third, an axis to plot the data as a colored grid
ax = fig.add_subplot(grid[3:22, 12:19])

# Define limits of colormap is none are provided, robust max and min using percentiles
if vmin is None and vmax is None:
    vmax = np.nanpercentile(df_sub[statistic_name].values, 99)
    vmin = np.nanpercentile(df_sub[statistic_name].values, 1)

# Create custom colormap
# col_bounds = np.array([vmin, np.mean([vmin, vmax]), vmax])
col_bounds = np.array([0.8, 1.2, 2, 3.5, 5])
cb = []
cb_val = np.linspace(0, 1, len(col_bounds))
for j in range(len(cb_val)):
    cb.append(cmap(cb_val[j]))
cmap_cus = colors.LinearSegmentedColormap.from_list('my_cb', list(
    zip((col_bounds - min(col_bounds)) / (max(col_bounds - min(col_bounds))), cb)), N=1000)

# Plot a 2D colored grid using fill_between
list_count_mismatch = []
for i in range(len(unique_var_1)):
    for j in range(len(unique_var_2)):
        df_both = df_sub[np.logical_and(df_sub.var1_mid == unique_var_1[i], df_sub.var2_mid == unique_var_2[j])]
        df_both_gla = df_sub_gla[np.logical_and(df_sub_gla.var1_mid == unique_var_1[i], df_sub_gla.var2_mid == unique_var_2[j])]
        df_both_all = df_sub_all[np.logical_and(df_sub_all.var1_mid == unique_var_1[i], df_sub_all.var2_mid == unique_var_2[j])]

        stat = df_both[statistic_name].values[0]
        stat_gla = df_both_gla[statistic_name].values[0]
        stat_all = df_both_all[statistic_name].values[0]

        diff_stable_gla = stat_gla - stat

        if np.isfinite(stat_all):
            stat_col_all = max(0.0001, min(0.9999, (stat_all - min(col_bounds)) / (max(col_bounds) - min(col_bounds))))
            col_all = cmap_cus(stat_col_all)
        else:
            col_all = nodata_color

        ax.fill_between([i-0.5, i+0.5], [j-0.5] * 2, [j+0.5] * 2, facecolor=col_all, alpha=1, edgecolor='white')

        if np.abs(diff_stable_gla/stat)<0.15:
            col_point = 'white'
        elif np.abs(diff_stable_gla/stat)<0.3:
            col_point='lightgrey'
        elif np.abs(diff_stable_gla / stat) < 0.5:
            col_point='darkgrey'
        else:
            col_point='black'

        if np.isfinite(stat_all):
            ax.add_patch(
                mpatches.Circle(xy=(i, j), radius=0.1, facecolor=col_point, edgecolor='white', alpha=1, zorder=30))

        if np.abs(diff_stable_gla/stat)>0.3:
            list_count_mismatch.append(df_both['count'].values)

        # ax.fill_between([i, i+0.5], [j] * 2, [j+0.5] * 2, facecolor=col, alpha=1, edgecolor='tab:brown')
        # ax.fill_between([i-0.5, i+0.5], [j-0.5] * 2, [j+0.5]*2, facecolor = 'None', edgecolor='black')

ax.set_xlabel(label_var_name_1)
ax.set_ylabel(label_var_name_2, labelpad=-2.5)
ax.set_xscale(scale_var_1)
ax.set_yscale(scale_var_2)
# In case the axis value does not agree with the scale (e.g., 0 for log scale)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ax.set_xlim((-0.5, len(bins_slope)-1.5))
    ax.set_ylim((-0.5, len(bins_curv)-1.5))

ax.set_xticks(np.arange(0, len(bins_slope)-1))
ax.set_xticklabels(list_bin_name_2_slope, rotation=90)
ax.set_yticks(np.arange(0, len(bins_curv)-1))
ax.set_yticklabels(list_bin_name_curv)

# 3.4/ Fourth and finally, add a colormap and nodata color to the legend
axcmap = fig.add_subplot(grid[:3, 19:22])

# Remove ticks, labels, frame
axcmap.set_xticks([])
axcmap.set_yticks([])
axcmap.spines['top'].set_visible(False)
axcmap.spines['left'].set_visible(False)
axcmap.spines['right'].set_visible(False)
axcmap.spines['bottom'].set_visible(False)

# Create an inset axis to manage the scale of the colormap
# cbaxes = axcmap.inset_axes([0.25, 0.75, 1, 0.2], label='cmap')
#
# # Create colormap object and plot
# norm = colors.Normalize(vmin=min(col_bounds), vmax=max(col_bounds))
# sm = plt.cm.ScalarMappable(cmap=cmap_cus, norm=norm)
# sm.set_array([])
# cb = plt.colorbar(sm, cax=cbaxes, orientation='horizontal', extend='both', shrink=0.8)
# cb.ax.tick_params(width=0.5, length=2)
# cb.set_label(label_statistic)

# Create an inset axis to manage the scale of the nodata legend
nodata = axcmap.inset_axes([0.05, 0, 1.25, 2.8], label='nodata')

# Plot a dot
nodata.add_patch(mpatches.Circle(xy=(0.05, 0.3), radius=0.03, facecolor='white', edgecolor='black', alpha=1, zorder=30, linewidth=0.25))
nodata.text(0.125, 0.3, '<15%', ha='left', va='center')
nodata.add_patch(mpatches.Circle(xy=(0.05, 0.2), radius=0.03, facecolor='lightgrey', edgecolor='black', alpha=1, zorder=30, linewidth=0.25))
nodata.text(0.125, 0.2, '15–30%', ha='left', va='center')
nodata.add_patch(mpatches.Circle(xy=(0.65, 0.3), radius=0.03, facecolor='darkgrey', edgecolor='black', alpha=1, zorder=30, linewidth=0.25))
nodata.text(0.725, 0.3, '30–50%', ha='left', va='center')
nodata.add_patch(mpatches.Circle(xy=(0.65, 0.2), radius=0.03, facecolor='black', edgecolor='black', alpha=1, zorder=30, linewidth=0.25))
nodata.text(0.725, 0.2, '>50%', ha='left', va='center')

nodata.text(0.5, 0.1, 'Dispersion\ndifference\nstable to\nmoving', ha='center', va='top')
nodata.set_xlim((0, 1))
nodata.set_ylim((0, 1))
nodata.set_xticks([])
nodata.set_yticks([])
nodata.spines['top'].set_visible(False)
nodata.spines['left'].set_visible(False)
nodata.spines['right'].set_visible(False)
nodata.spines['bottom'].set_visible(False)

# # Plot a nodata legend
# nodata.fill_between([0, 1], [0, 0], [1, 1], facecolor=nodata_color)
#
# nodata.text(0.5, -0.25, 'No data', ha='center', va='top')

# 4/ Panel d: plot the resulting error map for the Mont-Blanc case study

# First, compute the error function from N-D interpolation

fn = xdem.spatialstats.interp_nd_binning(df_sub, list_var_names=['slope', 'maxc'], statistic='nmad', min_count=30)
maxabsc[maxabsc>50] = 50
dh_err = fn((slope.data, maxabsc))

# 4.1/ Plot the error map

ax = fig.add_subplot(grid[:-2, 24:], projection=ccrs.UTM(32), label='Mont-Blanc')

# ax.text(0.6, 0.975, dates[i], transform=ax.transAxes, fontweight='bold', va='top', ha='right',bbox=dict(facecolor='white', alpha=1))

dh_err = pleia_ddem.copy(new_array=dh_err)
dh_err_full = dh_err.copy()

y_extent = pleia_ddem.bounds.top - pleia_ddem.bounds.bottom
x_extent = pleia_ddem.bounds.right - pleia_ddem.bounds.left
plt_extent = [
    pleia_ddem.bounds.left + 3/10*x_extent,
    pleia_ddem.bounds.right - 2/10*x_extent,
    pleia_ddem.bounds.bottom + 1/10*y_extent,
    pleia_ddem.bounds.top - 4.3/10*y_extent,
]
crop_ext = [plt_extent[0], plt_extent[2], plt_extent[1], plt_extent[3]]
dh_err.crop(crop_ext)

# ax.add_feature(cfeature.NaturalEarthFeature('physical', 'ocean', '50m', facecolor='gainsboro'))
# ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', facecolor='grey'))
# ax.outline_patch.set_edgecolor('black')

hs = gu.Raster(fn_hs)
hs.crop(crop_ext)
hs_arr = hs.data
# def stretch_hs(hs,stretch_factor=1.):
#
#     max_hs = 255
#     min_hs = 0
#
#     hs_s = (hs - (max_hs-min_hs)/2)*stretch_factor + (max_hs-min_hs)/2
#
#     return hs_s
#
# hs_arr = stretch_hs(hs_arr ,stretch_factor=0.9)
#
# xx, yy = hs.coords()
#
# ind_ll = xx <= plt_extent[0] + 1/4 * x_extent
# ind_ml = np.logical_and(xx > plt_extent[0] + 1/4 * x_extent, xx <= plt_extent[0] + 1/2*x_extent)
# ind_mr = np.logical_and(xx > plt_extent[0] + 1/2 * x_extent, xx <= plt_extent[0] + 3/4*x_extent)
# ind_rr = xx > plt_extent[0] + 3/4 * x_extent
# list_inds = [ind_ll, ind_ml, ind_mr, ind_rr]
# list_alpha = [1, 1, 1, 1]
#
# for i, ind in enumerate(list_inds):
#     hs_arr_tmp = np.copy(hs_arr)
#     hs_arr_tmp[0, ~ind] = 0
#
hs_arr[~np.isfinite(dh_err.data)] = 0
color1 = colors.to_rgba('black')
color2 = colors.to_rgba('white')
cmap_ll = colors.LinearSegmentedColormap.from_list('my_cmap_hs', [color1, color2], 256)
cmap_ll._init()
cmap_ll._lut[1:, -1] = 1
cmap_ll._lut[0:1, -1] = 0.0  # We make transparent the lowest hillshade value

ax.imshow(hs_arr[0, :, :], extent=plt_extent, transform=ccrs.UTM(32), cmap=cmap_ll,
      interpolation=None, zorder=2)

ax.add_patch(mpatches.Rectangle((crop_ext[2] - 2.5*x_extent/10, crop_ext[1] + y_extent*39/40),2000, 200,
                                edgecolor='black',facecolor='black',transform=ccrs.UTM(32),zorder=10,linewidth=0.5))
ax.add_patch(mpatches.Rectangle((crop_ext[2] - 2.5*x_extent/10 + 2000, crop_ext[1] + y_extent*39/40),2000, 200,
                                edgecolor='black',facecolor='white',transform=ccrs.UTM(32),zorder=10,linewidth=0.5))
ax.text(crop_ext[2] - 2.5*x_extent/10,  crop_ext[1] + y_extent*39/40 - 100,'0 km',ha='center',va='top',transform=ccrs.UTM(32),zorder=10)
ax.text(crop_ext[2] - 2.5*x_extent/10 + 2000, crop_ext[1] + y_extent*39/40 - 100,'2 km',ha='center',va='top',transform=ccrs.UTM(32),zorder=10)
ax.text(crop_ext[2] - 2.5*x_extent/10 + 4000, crop_ext[1] + y_extent*39/40 - 100,'4 km',ha='center',va='top',transform=ccrs.UTM(32),zorder=10)

# cmap2._lut[0, -1] = 0.0  # We made transparent de 10 first levels of hillshade,
cmap_cus.set_bad(color='None')
ax.imshow(dh_err.data[0, :, :], extent=plt_extent, transform=ccrs.UTM(32), vmin=0.8, vmax=5, cmap=cmap_cus,
          interpolation=None, zorder=3, alpha=0.85)

shape_feature = ShapelyFeature(Reader(fn_shp).geometries(), ccrs.PlateCarree(), edgecolor=plt.cm.Greys(0.75), alpha=1,
                               facecolor='None', linewidth=0.75, zorder=4)
ax.add_feature(shape_feature)
ax.gridlines(draw_labels={'top':'x', 'right':'y'}, dms=True, x_inline=False, y_inline=False)

ax.text(0.025, 0.975, 'd', transform=ax.transAxes, ha='left', va='top',fontweight='bold', fontsize=14, zorder=30)


# 4.2/ Add an inset that zooms on the area of Figure 1 panel a

inside_ext = [338680, 5086760, 340680, 5087460]
plt_extent = [inside_ext[0], inside_ext[2], inside_ext[1], inside_ext[3]]

ax.add_patch(mpatches.Rectangle((inside_ext[0], inside_ext[1]), inside_ext[2] - inside_ext[0], inside_ext[3] - inside_ext[1],
                                edgecolor='black',facecolor='None',transform=ccrs.UTM(32),
                                zorder=10, linewidth=1.5))

inside_ax = fig.add_axes([0.77, 0.05, 0.15, 0.15], projection=ccrs.UTM(32), label='subset')


inside_ax.spines['geo'].set_edgecolor('black')
inside_ax.spines['geo'].set_linewidth(1.5)
dh_err_full.crop(inside_ext)
inside_ax.imshow(dh_err_full.data[0, :, :], extent=plt_extent, transform=ccrs.UTM(32), vmin=0.8, vmax=5, cmap=cmap_cus,
          interpolation=None, zorder=3, alpha=0.85)

shape_feature = ShapelyFeature(Reader(fn_shp).geometries(), ccrs.PlateCarree(), edgecolor=plt.cm.Greys(0.75), alpha=1,
                               facecolor='None', linewidth=0.75, zorder=4)
inside_ax.add_feature(shape_feature)
# inside_ax.gridlines(draw_labels=False, dms=False, x_inline=False, y_inline=False)

hs = gu.Raster(fn_hs)
hs.crop(inside_ext)
hs_arr = hs.data

inside_ax.imshow(hs_arr[0, :, :], extent=[inside_ext[0], inside_ext[2], inside_ext[1], inside_ext[3]],
                 transform=ccrs.UTM(32), cmap=cmap_ll, interpolation=None, zorder=2)


ax.add_patch(mpatches.ConnectionPatch(xyA=(inside_ext[0], inside_ext[1]), xyB=(0, 1), coordsA=ax.transData, coordsB=inside_ax.transAxes,
                                      zorder=10, facecolor='black'))
ax.add_patch(mpatches.ConnectionPatch(xyA=(inside_ext[2], inside_ext[1]), xyB=(1, 1), coordsA=ax.transData, coordsB=inside_ax.transAxes,
                                      zorder=10, facecolor='black'))


# 4.3/ Add colormap and legend

cbaxes = ax.inset_axes([-0.33, -0.05, 0.04, 0.75], zorder=10)

norm = colors.Normalize(vmin=0.8, vmax=5)
sm = plt.cm.ScalarMappable(cmap=cmap_cus, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, cax=cbaxes, ticks=[0.8, 2, 5], orientation='vertical', extend='both', shrink=0.2)
cb.set_label('Dispersion (1$\sigma$) of elevation differences (m)\n inferred from stable terrain')


axleg = ax.inset_axes([0.1, -0.125, 0.1, 0.1], zorder=1)
axleg.set_xticks([])
axleg.set_yticks([])
axleg.spines['top'].set_visible(False)
axleg.spines['left'].set_visible(False)
axleg.spines['right'].set_visible(False)
axleg.spines['bottom'].set_visible(False)

axleg.add_patch(mpatches.Rectangle((0.2, 0.2), 0.7, 0.7, edgecolor=plt.cm.Greys(0.75),
                                linewidth=0.75, facecolor='None', alpha=1))
axleg.set_xlim((0, 1))
axleg.set_ylim((0, 1))
axleg.text(0.55, 0.02, 'Glacier\noutlines', ha='center', va='top')

# Save to file
plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_4_final.png', dpi=400)

# Additionally, compute the percentage of terrain with more than 10% difference between the variance of stable and unstable
s0=np.sum(list_count_mismatch)
df_tmp = df_sub_all[df_sub_all.nd==2]
s1=np.sum(df_tmp['count'])
perc = s0/s1