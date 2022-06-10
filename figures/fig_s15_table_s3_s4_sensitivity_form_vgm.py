"""Plotting of Figure S15 and Table S3-4: sensitivity with the form of variogram models for the Mont-Blanc case study"""
import numpy as np
import pandas as pd
import xdem
import matplotlib.pyplot as plt

# Open empirical variograms
df_std_sta = pd.read_csv('/home/atom/ongoing/work_stderr_dem/case_study_montblanc/df_vgm_std_sta.csv')

df_all_patches = pd.read_csv('/home/atom/ongoing/work_stderr_dem/case_study_montblanc/df_patches_sta_gla.csv')

df_std_sta = df_std_sta[df_std_sta.bins<30000]
df_std_sta.err_exp /= np.sqrt(100)/2

std_fac_sta = np.nanmean(df_std_sta.exp.values[-3:])
df_std_sta.exp /= std_fac_sta
df_std_sta.err_exp /= std_fac_sta

name_types = ['spherical', 'exponential', 'gaussian']
list_types = ['Sph', 'Exp', 'Gau']
list_fun, list_params, list_labels = ([] for i in range(3))
for i in range(len(list_types)):
    for j in range(len(list_types)):

        fun, params = xdem.spatialstats.fit_sum_model_variogram(list_model=[list_types[i], list_types[j]], empirical_variogram=df_std_sta,
                                                             bounds=[(0, 1000), (0, 1), (1000, 20000), (0, 1)],
                                                             p0=[500, 0.5, 15000, 0.5])

        list_fun.append(fun)
        list_params.append(params)
        list_labels.append('Short: '+name_types[i]+'\nLong: '+name_types[j])


# Initiate figure
fig = plt.figure(figsize=(12,5.5))

ylabel = "Standardized\nelevation variance"
xlabel = 'Spatial lag (m)'
xscale_range_split = [120, 640, 3600]
xscale = 'linear'
list_fit_fun = list_fun
list_fit_fun_label = list_labels
xlim = None
ylim = (0, 1.15)
list_df = [df_std_sta]
col_df = ['tab:orange', 'tab:red', 'tab:brown', 'tab:blue', 'tab:green', 'tab:olive', 'lightgrey', 'darkgrey', 'black']
label_df = ['Stable', 'Unstable']
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
        ax0.fill_between([mid-width/3, mid+width/3], [0] * 2, [count] * 2,
                         facecolor='black', alpha=1,
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
        ax.errorbar(bins_center, df.exp, yerr=df.err_exp, fmt='x', color='black')

    for i in range(len(df0)):
        width = interval_var[i + 1] - interval_var[i]
        mid = interval_var[i] + width / 2
        ax.vlines(mid - width / 2, ymin=[0], ymax=2*max(df0.exp), colors='tab:gray', linestyles='dashed', linewidths=0.5)

    # ax.hlines(1, xmin=xmin[k], xmax=xmax[k], colors='black', linestyles='dotted')
    # If a list of functions is passed, plot the modelled variograms
    if list_fit_fun is not None:
        for i, fit_fun in enumerate(list_fit_fun):
            x = np.linspace(xmin[k], xmax[k], 1000)

            ax.plot(x, fit_fun(x), linestyle='dashdot', color=col_df[i], zorder=30, linewidth=1.5, label=list_fit_fun_label[i])

    ax.set_xscale(xscale)
    ax.set_xticks([])

    if xlim is None:
        ax.set_xlim((xmin[k], xmax[k]))
    else:
        ax.set_xlim(xlim)

    if ylim is not None:
        ax.set_ylim(ylim)
    else:
        ax.set_ylim((0, np.nanmax(df_std_sta.exp) + np.nanmean(df_std_sta.err_exp)))

    if k == nb_subpanels - 1:
        handles, labels = plt.gca().get_legend_handles_labels()
        order = [0, 1, 2, 3, 4, 5, 6, 7, 8]
        ax.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='lower right', ncol=3)
    if k == 0:
        ax.set_ylabel(ylabel)
        ax.text(0.1, 0.95, 'a', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)

    else:
        ax.set_yticks([])


ylim = (0.9, 1.025)
# Loop over each subpanel
for k in [1,0,2,3]:
    # Now, plot the statistic of the data
    ax = fig.add_subplot(grid[11:, slice(xgridmin[k],xgridmax[k])])

    # Get the bins center
    for i, df in enumerate(list_df):
        bins_center = np.subtract(df.bins, np.diff([0] + df.bins.tolist()) / 2)
        ax.scatter(bins_center, df.exp, marker='x', color='black')

    for i in range(len(df0)):
        width = interval_var[i + 1] - interval_var[i]
        mid = interval_var[i] + width / 2
        ax.vlines(mid - width / 2, ymin=[0], ymax=2*max(df0.exp), colors='tab:gray', linestyles='dashed', linewidths=0.5)

    # ax.hlines(1, xmin=xmin[k], xmax=xmax[k], colors='black', linestyles='dotted')
    # If a list of functions is passed, plot the modelled variograms
    if list_fit_fun is not None:
        for i, fit_fun in enumerate(list_fit_fun):
            x = np.linspace(xmin[k], xmax[k], 1000)

            ax.plot(x, fit_fun(x), linestyle='dashdot', color=col_df[i], zorder=30, linewidth=1.5, label=list_fit_fun_label[i])

    ax.set_xscale(xscale)
    if k == 1:
        ax.set_xlabel(xlabel, x=1, ha='center')

    # ax.set_xticks([])

    if xlim is None:
        ax.set_xlim((xmin[k], xmax[k]))
    else:
        ax.set_xlim(xlim)

    if ylim is not None:
        ax.set_ylim(ylim)
    else:
        ax.set_ylim((0, np.nanmax(df_std_sta.exp) + np.nanmean(df_std_sta.err_exp)))

    if k == nb_subpanels - 1:
        ax.scatter([], [], color='black', label='Empirical variogram', marker='x')
        handles, labels = plt.gca().get_legend_handles_labels()
        order = [9]
        ax.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='lower right', ncol=1)
    if k == 0:
        ax.set_ylabel(ylabel)
        ax.text(0.1, 0.95, 'b', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)

    else:
        ax.set_yticks([])


# Save to file
plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_S15_final.png', dpi=400)

# Tables S3 and S4
list_rss = [sum((f(df_std_sta.bins.values) - df_std_sta.exp.values)**2) for f in list_fun]
list_rss_long = [sum((f(df_std_sta.bins.values[df_std_sta.bins.values>1000]) - df_std_sta.exp.values[df_std_sta.bins.values>1000])**2) for f in list_fun]
list_rss_short = [sum((f(df_std_sta.bins.values[df_std_sta.bins.values<=1000]) - df_std_sta.exp.values[df_std_sta.bins.values<1000])**2) for f in list_fun]