"""Plotting of Figure S6: heteroscedasticity with slope and quality of stereo-correlation for the NPI case study"""
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
import matplotlib.colors as colors
from matplotlib.legend_handler import HandlerLineCollection, HandlerPatch

# Same script as for main Figure 4, adapted for the NPI case study
fn_ddem = '/home/atom/ongoing/work_stderr_dem/case_study_npi/dh_ASTER-SPOT_NK_Deramp.tif'
fn_shp = '/home/atom/data/inventory_products/RGI/00_rgi60_neighb_merged/17_rgi60_SouthernAndes/17_rgi60_SouthernAndes.shp'
fn_corr='/home/atom/ongoing/work_stderr_dem/case_study_npi/AST_L1A_00303182012144228/AST_L1A_00303182012144228_CORR.tif'
fn_dem_aster = '/home/atom/ongoing/work_stderr_dem/case_study_npi/AST_L1A_00303182012144228/AST_L1A_00303182012144228_Z.tif'

pleia_ddem = gu.Raster(fn_ddem)
ref_dem = gu.Raster(fn_dem_aster)
glaciers_outlines = gu.Vector(fn_shp)
mask_glacier = glaciers_outlines.create_mask(pleia_ddem)

# Remove forest, very large outliers
pleia_ddem.data[np.abs(pleia_ddem.data)>200] = np.nan

slope, planc, profc = xdem.terrain.get_terrain_attribute(ref_dem, attribute=['slope', 'planform_curvature',
                                                                                'profile_curvature'])

corr = gu.Raster(fn_corr)
maxabsc, _ = gu.spatial_tools.get_array_and_mask(corr)
maxabsc = maxabsc[None, :, :]

# Filter large outliers per category
bins_slope = [0, 2.5, 5, 10, 15, 20, 30, 40, 50, 70, 90]
# bins_curv = [-10, -8, -6, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 6, 8, 10]
# bins_curv = [0, 0.2, 0.5, 1, 2, 3, 4, 6, 10, 20, 50]
bins_curv = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

for i in range(len(bins_slope) - 1):
    # Subset by slope category
    subset = np.logical_and(slope.data >= bins_slope[i], slope.data < bins_slope[i + 1])
    dh_sub = pleia_ddem.data[subset]
    # Remove very large outliers of the category
    med_sub = np.nanmedian(dh_sub)
    nmad_sub = xdem.spatialstats.nmad(dh_sub)
    pleia_ddem.data[np.logical_and(subset, np.abs(pleia_ddem.data-med_sub) > 5 * nmad_sub)] = np.nan
for i in range(len(bins_curv) - 1):
    # Subset by slope category
    subset = np.logical_and(maxabsc >= bins_curv[i], maxabsc < bins_curv[i + 1])
    dh_sub = pleia_ddem.data[subset]
    # Remove very large outliers of the category
    med_sub = np.nanmedian(dh_sub)
    nmad_sub = xdem.spatialstats.nmad(dh_sub)
    pleia_ddem.data[np.logical_and(subset, np.abs(pleia_ddem.data-med_sub) > 5 * nmad_sub)] = np.nan

df = xdem.spatialstats.nd_binning(pleia_ddem.data.ravel(), list_var=[slope.data.ravel()], list_var_names=['slope'])

# Format data for seaborn plot
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

# Create custom colormap
# col_bounds = np.array([vmin, np.mean([vmin, vmax]), vmax])
cmap = plt.get_cmap('YlOrRd')
col_bounds = np.array([5, 10, 20, 30, 50])
cb = []
cb_val = np.linspace(0, 1, len(col_bounds))
for j in range(len(cb_val)):
    cb.append(cmap(cb_val[j]))
cmap_cus = colors.LinearSegmentedColormap.from_list('my_cb', list(
    zip((col_bounds - min(col_bounds)) / (max(col_bounds - min(col_bounds))), cb)), N=1000)

# Splitted violin plot
fig = plt.figure(figsize=(15, 7))
grid = plt.GridSpec(24, 24, wspace=0.1, hspace=0.1)

# First, an horizontal axis on top to plot the sample histograms
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

# Then, another axis to plot the splitted violin plots
ax = fig.add_subplot(grid[2:9, :10])


# Plot a colored curve for the dispersion

# Get a linearized version of the bins to convert from discrete to continuous plotting
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
lc = LineCollection(segments, cmap=cmap_cus, zorder=2)
lc.set_array(yfit)
lc.set_linewidth(5)
ax.add_collection(lc)

a = list(ax.get_children())

sns.violinplot(ax=ax, x="slope_category", y="dh", hue="terrain",
               data=df_all, split=True,  palette={'Moving':'tab:cyan', 'Stable':'tab:brown'}, linewidth=1,
               edgecolor='black')

ax.set_ylim((-95, 95))

del df_all

b = list(ax.get_children())
c = set(b) - set(a)
for d in c:
    d.set_zorder(1000)

ax.set_xlabel('Slope categories (degrees)')
ax.set_ylabel('Elevation\ndifferences\n(m)', labelpad=-5)
ax.set_xlim((-0.5, len(list_counts)-0.5))


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

# Format data for seaborn plot
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
    bin_name =  '{:.0f}'.format(bins_curv[i])+'–\n'+'{:.0f}'.format(bins_curv[i+1])


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

# First, an horizontal axis on top to plot the sample histograms
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

# Then, another axis to plot the splitted violin plots
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

ax.set_ylim((-95, 95))

del df_all_curv

b = list(ax.get_children())
c = set(b)-set(a)
for d in c:
    d.set_zorder(1000)

ax.set_xlabel('Quality of stereo-correlation (%)')
ax.set_ylabel('Elevation\ndifferences\n(m)',labelpad=-5)
leg = ax.legend(loc='lower left')
leg.remove()
ax.set_xlim((-0.5, len(list_counts)-0.5))

## Finally, plot the carpet with double distribution
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
label_var_name_2 = 'Quality of stereo-correlation (%)'

label_statistic = 'Dispersion (1$\sigma$) of\nelevation\ndifferences (m)'

nodata_color = 'white'
from matplotlib import colors

# First, an horizontal axis on top to plot the sample histogram of the first variable
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


# Second, a vertical axis on the right to plot the sample histogram of the second variable
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

# Third, an axis to plot the data as a colored grid
ax = fig.add_subplot(grid[3:22, 12:19])

# Define limits of colormap is none are provided, robust max and min using percentiles
if vmin is None and vmax is None:
    vmax = np.nanpercentile(df_sub[statistic_name].values, 99)
    vmin = np.nanpercentile(df_sub[statistic_name].values, 1)

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

        if np.abs(diff_stable_gla/stat)>0.30:
            list_count_mismatch.append(df_both['count'].values)

        # ax.fill_between([i, i+0.5], [j] * 2, [j+0.5] * 2, facecolor=col, alpha=1, edgecolor='tab:brown')
        # ax.fill_between([i-0.5, i+0.5], [j-0.5] * 2, [j+0.5]*2, facecolor = 'None', edgecolor='black')

ax.set_xlabel(label_var_name_1)
ax.set_ylabel(label_var_name_2)
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


cbaxes = ax.inset_axes([1.25, 0.2, 0.04, 0.6], zorder=10)

norm = colors.Normalize(vmin=5, vmax=50)
sm = plt.cm.ScalarMappable(cmap=cmap_cus, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, cax=cbaxes, ticks=[5, 10, 20, 50], orientation='vertical', extend='both', shrink=0.2)
cb.set_label('Dispersion (1$\sigma$) of elevation differences (m)\n inferred from stable terrain')


# Fourth and finally, add a colormap and nodata color to the legend
axcmap = fig.add_subplot(grid[:3, 19:22])

# Remove ticks, labels, frame
axcmap.set_xticks([])
axcmap.set_yticks([])
axcmap.spines['top'].set_visible(False)
axcmap.spines['left'].set_visible(False)
axcmap.spines['right'].set_visible(False)
axcmap.spines['bottom'].set_visible(False)

# Create an inset axis to manage the scale of the nodata legend
nodata = axcmap.inset_axes([0.05, 0.1, 1.25, 2.8], label='nodata')

# Plot a dot
nodata.add_patch(mpatches.Circle(xy=(0.05, 0.3), radius=0.03, facecolor='white', edgecolor='black', alpha=1, zorder=30, linewidth=0.25))
nodata.text(0.125, 0.3, '<15%', ha='left', va='center')
nodata.add_patch(mpatches.Circle(xy=(0.05, 0.2), radius=0.03, facecolor='lightgrey', edgecolor='black', alpha=1, zorder=30, linewidth=0.25))
nodata.text(0.125, 0.2, '15–30%', ha='left', va='center')
nodata.add_patch(mpatches.Circle(xy=(0.65, 0.3), radius=0.03, facecolor='darkgrey', edgecolor='black', alpha=1, zorder=30, linewidth=0.25))
nodata.text(0.725, 0.3, '30–50%', ha='left', va='center')
nodata.add_patch(mpatches.Circle(xy=(0.65, 0.2), radius=0.03, facecolor='black', edgecolor='black', alpha=1, zorder=30, linewidth=0.25))
nodata.text(0.725, 0.2, '>50%', ha='left', va='center')

nodata.text(0.5, 0.1, 'Dispersion difference\nstable to moving', ha='center', va='top')
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

plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_S6_final.png', dpi=300)

# Percentage of terrain with more than 10% difference between the variance of stable and unstable
s0=np.sum(list_count_mismatch)
df_tmp = df_sub_all[df_sub_all.nd==2]
s1=np.sum(df_tmp['count'])
perc = s0/s1