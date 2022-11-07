"""Plotting of Figure S11: effect of variogram estimator robustness for the Mont-Blanc case study"""
import numpy as np
import pandas as pd
import xdem
import geoutils as gu
import matplotlib.pyplot as plt
import skgstat as skg

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
pleia_ddem.data[np.abs(pleia_ddem.data)>500] = np.nan

pleia_ddem.data[mask_glacier] = np.nan

slope, planc, profc = xdem.terrain.get_terrain_attribute(ref_dem, attribute=['slope', 'planform_curvature',
                                                                                'profile_curvature'])
maxabsc = np.maximum(np.abs(planc), np.abs(profc))

# First, without outlier filtering

fn_dowd_nofilt = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/df_vgm_sta_dowd_nofilt.csv'
fn_matheron_nofilt = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/df_vgm_sta_matheron_nofilt.csv'

# df_vgm_sta_dowd_nofilt = xdem.spatialstats.sample_empirical_variogram(pleia_ddem.data.data, pleia_ddem.res[0], subsample=50, n_variograms=10, runs=20, estimator='dowd', n_jobs=3)
# df_vgm_sta_matheron_nofilt = xdem.spatialstats.sample_empirical_variogram(pleia_ddem.data.data, pleia_ddem.res[0], subsample=50, n_variograms=10, runs=20, estimator='matheron', n_jobs=3)
# df_vgm_sta_dowd_nofilt.to_csv(fn_dowd_nofilt, index=False)
# df_vgm_sta_matheron_nofilt.to_csv(fn_matheron_nofilt, index=False)

df_vgm_sta_dowd_nofilt = pd.read_csv(fn_dowd_nofilt)
# Normalizing to semi-variance
df_vgm_sta_dowd_nofilt.exp /= 2
df_vgm_sta_dowd_nofilt.err_exp /= 2
df_vgm_sta_matheron_nofilt = pd.read_csv(fn_matheron_nofilt)

# # Filter large outliers per category
bins_slope = [0, 2.5, 5, 10, 15, 20, 30, 40, 50, 70, 90]
# # bins_curv = [-10, -8, -6, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 6, 8, 10]
bins_curv = [0, 0.2, 0.5, 1, 2, 3, 4, 6, 10, 20, 50]
for i in range(len(bins_slope) - 1):
    # Subset by slope category
    subset = np.logical_and(slope.data >= bins_slope[i], slope.data < bins_slope[i + 1])
    dh_sub = pleia_ddem.data[subset]
    # Remove very large outliers of the category
    med_sub = np.nanmedian(dh_sub)
    nmad_sub = xdem.spatialstats.nmad(dh_sub)
    # Remove outliers
    pleia_ddem.data[np.logical_and(subset, np.abs(pleia_ddem.data-med_sub) > 7 * nmad_sub)] = np.nan
for i in range(len(bins_curv) - 1):
    # Subset by slope category
    subset = np.logical_and(maxabsc >= bins_curv[i], maxabsc < bins_curv[i + 1])
    dh_sub = pleia_ddem.data[subset]
    # Remove very large outliers of the category
    med_sub = np.nanmedian(dh_sub)
    nmad_sub = xdem.spatialstats.nmad(dh_sub)
    # Remove outliers
    pleia_ddem.data[np.logical_and(subset, np.abs(pleia_ddem.data-med_sub) > 7 * nmad_sub)] = np.nan


# Then, with outlier filtering

fn_dowd_filt = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/df_vgm_sta_dowd_filt.csv'
fn_matheron_filt = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/df_vgm_sta_matheron_filt.csv'

# df_vgm_sta_dowd = xdem.spatialstats.sample_empirical_variogram(pleia_ddem.data.data, pleia_ddem.res[0], subsample=50, n_variograms=10, runs=20, estimator='dowd', n_jobs=3)
# df_vgm_sta_matheron = xdem.spatialstats.sample_empirical_variogram(pleia_ddem.data.data, pleia_ddem.res[0], subsample=50, n_variograms=10, runs=20, estimator='matheron', n_jobs=3)
# df_vgm_sta_dowd.to_csv(fn_dowd_filt, index=False)
# df_vgm_sta_matheron.to_csv(fn_matheron_filt, index=False)

df_vgm_sta_dowd = pd.read_csv(fn_dowd_filt)
# Normalizing to semi-variance
df_vgm_sta_dowd.exp /= 2
df_vgm_sta_dowd.err_exp /= 2
df_vgm_sta_matheron = pd.read_csv(fn_matheron_filt)


# FIGURE
fig = plt.figure(figsize=(12,5.5))

ylabel = "Elevation variance (m²)"
xlabel = 'Spatial lag (m)'
xscale_range_split = [120, 640, 3600]
xscale = 'linear'
list_fit_fun = None
list_fit_fun_label = None
xlim = None
ylim = (0, 55)
list_df = [df_vgm_sta_matheron_nofilt, df_vgm_sta_dowd_nofilt]
col_df = ['tab:orange', 'tab:blue']
label_df = ['Matheron estimator', 'Dowd estimator']
df0 = list_df[0]

init_gridsize = [10, 20]
# Create parameters to split x axis into different linear scales
# If there is no split, get parameters for a single subplot
if xscale_range_split is None:
    nb_subpanels = 1
    if xscale == 'log':
        xmin = [np.min(df0.bins) / 2]
    else:
        xmin = [0]
    xmax = [np.max(df0.bins)]
    xgridmin = [0]
    xgridmax = [init_gridsize[0]]
    gridsize = init_gridsize
# Otherwise, derive a list for each subplot
else:
    # Add initial zero if not in input
    if xscale_range_split[0] != 0:
        if xscale == 'log':
            first_xmin = np.min(df0.bins) / 2
        else:
            first_xmin = 0
        xscale_range_split = [first_xmin] + xscale_range_split
    # Add maximum distance if not in input
    if xscale_range_split[-1] != np.max(df0.bins):
        xscale_range_split.append(15000)

    # Scale grid size by the number of subpanels
    nb_subpanels = len(xscale_range_split) - 1
    gridsize = init_gridsize.copy()
    gridsize[1] = 18
    gridsize[0] *= nb_subpanels
    # Create list of parameters to pass to ax/grid objects of subpanels
    xmin, xmax, xgridmin, xgridmax = ([] for i in range(4))
    for i in range(nb_subpanels):
        xmin.append(xscale_range_split[i])
        xmax.append(xscale_range_split[i + 1])
        xgridmin.append(init_gridsize[0] * i)
        xgridmax.append(init_gridsize[0] * (i + 1))

# Need a grid plot to show the sample count and the statistic
grid = plt.GridSpec(gridsize[1], gridsize[0], wspace=0.1, hspace=0.1)

# Loop over each subpanel
for k in [1,0,2,3]:
    # First, an axis to plot the sample histogram
    ax0 = fig.add_subplot(grid[:3, xgridmin[k]:xgridmax[k]])
    ax0.set_xscale(xscale)
    ax0.set_xticks([])

    # Plot the histogram manually with fill_between
    interval_var = [0] + list(df0.bins)
    for i in range(len(df0)):
        width = interval_var[i+1] - interval_var[i]
        mid = interval_var[i] + width/2
        count = list_df[0]['count'].values[i]
        count_gla = list_df[1]['count'].values[i]
        ax0.fill_between([mid-width/3, mid], [0] * 2, [count] * 2,
                         facecolor='tab:orange', alpha=1,
                         edgecolor='black', linewidth=0.5)
        ax0.fill_between([mid, mid+width/3], [0] * 2, [count_gla] * 2,
                         facecolor='tab:blue', alpha=1,
                         edgecolor='black', linewidth=0.5)
        ax0.vlines(mid-width/2, ymin=[0], ymax=1.2*max(list_df[0]['count'].values), colors='tab:gray', linestyles='dashed', linewidths=0.5)
    if k == 0:
        ax0.set_ylabel('Pairwise\nsample\ncount')
        # Scientific format to avoid undesired additional space on the label side
        ax0.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    else:
        ax0.set_yticks([])
    # Ignore warnings for log scales
    ax0.set_xlim((xmin[k], xmax[k]))
    ax0.set_ylim((0, 1.2*max(list_df[0]['count'].values)))

    # Now, plot the statistic of the data
    ax = fig.add_subplot(grid[3:10, slice(xgridmin[k],xgridmax[k])])

    # Get the bins center
    for i, df in enumerate(list_df):
        bins_center = np.subtract(df.bins, np.diff([0] + df.bins.tolist()) / 2)
        ax.errorbar(bins_center, df.exp, yerr=df.err_exp, fmt='x', color=col_df[i])

    for i in range(len(df0)):
        width = interval_var[i + 1] - interval_var[i]
        mid = interval_var[i] + width / 2
        ax.vlines(mid - width / 2, ymin=[0], ymax=2*max(df0.exp), colors='tab:gray', linestyles='dashed', linewidths=0.5)

    # ax.hlines(1, xmin=xmin[k], xmax=xmax[k], colors='black', linestyles='dotted')
    # If a list of functions is passed, plot the modelled variograms
    if list_fit_fun is not None:
        for i, fit_fun in enumerate(list_fit_fun):
            x = np.linspace(xmin[k], xmax[k], 1000)

            def vgm_short(h):
                fn = skg.models.spherical(h, fit_fun[0], fit_fun[1])
                return fn

            def vgm_long(h):
                fn = skg.models.spherical(h, fit_fun[2], fit_fun[3])
                return fn

            def vgm_sum(h):
                fn = skg.models.spherical(h, fit_fun[0], fit_fun[1]) + skg.models.spherical(h, fit_fun[2], fit_fun[3])
                return fn

            colors_terrain = ['tab:brown', 'tab:cyan']

            ax.plot(x, vgm_short(x), linestyle='dashdot', color=colors_terrain[i], zorder=30, linewidth=1)
            ax.plot(x, vgm_sum(x), linestyle='dashed', color=colors_terrain[i], zorder=30, linewidth=1.5)
            # ax.plot(x, vgm_long(x), linestyle='dashdot', color=colors[i], label = 'Long-range model', zorder=30, linewidth=1)
            if i == 0:
                ax.errorbar([], [], [], color='black', label='Empirical variogram', fmt='x')
                ax.plot([], [], linestyle='dashdot', color='black', label='Modelled variogram: short-range')
                ax.plot([], [], linestyle='dashed', color='black', label='Modelled variogram: short- and long-range',
                        linewidth=1.5)
                ax.plot([], [], color='tab:brown', label='Stable terrain')
                ax.plot([], [], color='tab:cyan', label='Unstable terrain')


    ax.set_xscale(xscale)
    ax.set_xticks([])

    if xlim is None:
        ax.set_xlim((xmin[k], xmax[k]))
    else:
        ax.set_xlim(xlim)

    if ylim is not None:
        ax.set_ylim(ylim)
    else:
        ax.set_ylim((0, np.nanmax(df_vgm_sta_matheron_nofilt.exp) + np.nanmean(df_vgm_sta_matheron_nofilt.err_exp)))

    if k == nb_subpanels - 1:
        ax.errorbar([], [], [], color='tab:orange', label="Matheron's estimator with raw data", fmt='o')
        ax.errorbar([], [], [], color='tab:blue',
                    label="Dowd's estimator with raw data", fmt='o')
        # handles, labels = plt.gca().get_legend_handles_labels()
        # order = [4, 0, 1, 2, 3]
        # ax.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='lower right', ncol=2)
        ax.legend(loc='center right')
    if k == 0:
        ax.set_ylabel(ylabel)
        ax.text(0.1, 0.95, 'a', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)

    else:
        ax.set_yticks([])


list_df = [df_vgm_sta_matheron, df_vgm_sta_dowd]
list_fit_fun = None
col_df = ['tab:orange', 'tab:blue']
label_df = ['Matheron estimator', 'Dowd estimator']
ylim = (0, 11)


for k in [2,3,0,1]:

    # Now, plot the statistic of the data
    ax = fig.add_subplot(grid[10:, slice(xgridmin[k],xgridmax[k])])

    # Get the bins center
    list_fmt = ['o', 'o']
    for i, df in enumerate(list_df):
        if i == 0:
            bins_center = np.subtract(df.bins, np.diff([0] + df.bins.tolist()) / 2)
        else:
            bins_center = np.subtract(df.bins, np.diff([0] + df.bins.tolist()) / 2) + (xmax[k]-xmin[k])/20
        ax.errorbar(bins_center, df.exp, yerr=df.err_exp, fmt=list_fmt[i], color=col_df[i])

    for i in range(len(df0)):
        width = interval_var[i + 1] - interval_var[i]
        mid = interval_var[i] + width / 2
        ax.vlines(mid - width / 2, ymin=[0], ymax=2*max(df0.exp), colors='tab:gray', linestyles='dashed', linewidths=0.5)

    # ax.hlines(1, xmin=xmin[k], xmax=xmax[k], colors='black', linestyles='dotted')
    # If a list of functions is passed, plot the modelled variograms
    list_linestyle = ['dashed', 'dashdot']
    if list_fit_fun is not None:
        for i, fit_fun in enumerate(list_fit_fun):
            x = np.linspace(xmin[k], xmax[k], 1000)

            def vgm_short(h):
                fn = skg.models.spherical(h, fit_fun[0], fit_fun[1])
                return fn

            def vgm_long(h):
                fn = skg.models.spherical(h, fit_fun[2], fit_fun[3])
                return fn

            def vgm_sum(h):
                fn = skg.models.spherical(h, fit_fun[0], fit_fun[1]) + skg.models.spherical(h, fit_fun[2], fit_fun[3])
                return fn

            colors_terrain = ['tab:brown', 'tab:brown']

            # ax.plot(x, vgm_short(x), linestyle='dashdot', color=colors_terrain[i], zorder=30, linewidth=1)
            ax.plot(x, vgm_sum(x), linestyle=list_linestyle[i], color=colors_terrain[i], zorder=30, linewidth=1.5)
            # ax.plot(x, vgm_long(x), linestyle='dashdot', color=colors[i], label = 'Long-range model', zorder=30, linewidth=1)
            if i == 0:
                ax.plot([], [], linestyle='dashed', color='black', label='Modelled variogram sampled\nfrom standardized $dh$',
                        linewidth=1.5)
                ax.plot([], [], linestyle='dashdot', color='black', label='Modelled variogram divided by\naverage variance sampled from $dh$',
                        linewidth=1.5)


    ax.set_xscale(xscale)
    if nb_subpanels > 1 and k == (nb_subpanels - 1):
        ax.xaxis.set_ticks(np.linspace(xmin[k], xmax[k], 3))
    elif nb_subpanels > 1:
        ax.xaxis.set_ticks(np.linspace(xmin[k], xmax[k], 3)[:-1])

    if xlim is None:
        ax.set_xlim((xmin[k], xmax[k]))
    else:
        ax.set_xlim(xlim)

    if ylim is not None:
        ax.set_ylim(ylim)
    else:
        ax.set_ylim((0, np.nanmax(df_vgm_sta_matheron.exp) + np.nanmean(df_vgm_sta_matheron.err_exp)))

    if k == nb_subpanels - 1:
        pass
        # ax.annotate(text='', xy=(12500, 0.93), xytext=(7000, 0.82), arrowprops=dict(arrowstyle='-|>', connectionstyle=None,
        #                      shrinkA=0, shrinkB=0,
        #                      patchA=None, patchB=None, linewidth=1, facecolor='black'), zorder=30)
        # ax.text(7000, 0.8, '30-50% smaller error\nof variogram estimation\nusing standardized $dh$,\nmostly at long ranges', ha='center', va='top',
        #         bbox= dict(facecolor='white', boxstyle='round', alpha=0.8, linewidth=0.5), zorder=31)
    if k == 3:
        ax.errorbar([], [], [], color='tab:orange', label="Matheron's estimator with filtered data", fmt='o')
        ax.errorbar([], [], [], color='tab:blue',
                    label="Dowd's estimator with filtered data", fmt='o')

        # handles, labels = plt.gca().get_legend_handles_labels()
        # order = [0, 1, 2, 3]
        # ax.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='lower left', bbox_to_anchor=(-0.75, 0.05), ncol=2)
        ax.legend(loc='center right')
    if k == 1:
        ax.set_xlabel('Spatial lag (m)', x=1, ha='center')
    if k == 0:
        ax.set_ylabel(ylabel)
        ax.text(0.1, 0.95, 'b', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)
    else:
        ax.set_yticks([])

plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_S11_final.png', dpi=400)

