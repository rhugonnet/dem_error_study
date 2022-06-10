"""Plotting of Figure S8: heteroscedasticity interpolated from 2D slope/curvature binning for the Mont-Blanc case study"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import xdem

# Open file with estimates
fn_bin = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/df_heteroscedas_slope_curv.csv'

df_bin = pd.read_csv(fn_bin)

# Interpolate
fn = xdem.spatialstats.interp_nd_binning(df_bin, list_var_names=['slope_mid', 'maxc_mid'], statistic='nmad', min_count=30)

vec_slope = np.linspace(0, 90, 100)
vec_maxc = np.linspace(0, 15, 100)

grid = plt.GridSpec(22, 22, wspace=0.1, hspace=0.1)

# First, plot in linear scale
fig = plt.figure(figsize=(7,4))

ax = fig.add_subplot(grid[:, 0:9])

cmap = plt.get_cmap('YlOrRd')
col_bounds = np.array([0.8, 1.2, 2, 3.5, 5])
cb = []
cb_val = np.linspace(0, 1, len(col_bounds))
for j in range(len(cb_val)):
    cb.append(cmap(cb_val[j]))
cmap_cus = colors.LinearSegmentedColormap.from_list('my_cb', list(
    zip((col_bounds - min(col_bounds)) / (max(col_bounds - min(col_bounds))), cb)), N=1000)

for i in range(len(vec_slope)-1):
    for j in range(len(vec_maxc)-1):
        stat = fn([0.5*(vec_slope[i]+vec_slope[i+1]), 0.5*(vec_maxc[j]+vec_maxc[j+1])])
        if np.isfinite(stat):
            stat_col_all = max(0.0001, min(0.9999, (stat - min(col_bounds)) / (max(col_bounds) - min(col_bounds))))
            col = cmap_cus(stat_col_all)
        else:
            col = 'tab:gray'

        ax.fill_between(x=[vec_slope[i], vec_slope[i+1]], y1=[vec_maxc[j], vec_maxc[j]],
                        y2=[vec_maxc[j+1], vec_maxc[j+1]], facecolor=col)

ax.set_ylim((0, 15))
ax.set_xlim((0, 90))
ax.set_ylabel('Maximum absolute curvature (10$^{2}$ m$^{-1}$)')
ax.set_xlabel('Slope (degrees)')


# Create an inset axis to manage the scale of the colormap
cbaxes = ax.inset_axes([1.1, 0.2, 0.05, 0.6], label='cmap')

# Create colormap object and plot
norm = colors.Normalize(vmin=min(col_bounds), vmax=max(col_bounds))
sm = plt.cm.ScalarMappable(cmap=cmap_cus, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, cax=cbaxes, orientation='vertical', extend='both', shrink=0.8)
cb.ax.tick_params(width=0.5, length=2)
cb.set_label('Dispersion of elevation differences (1$\sigma$)')

ax.text(-0.025, 1.025, 'a', transform=ax.transAxes, ha='right', va='bottom',fontweight='bold', fontsize=14)


# Then, plot in log scales
ax = fig.add_subplot(grid[:, 13:])


vec_slope = np.exp(np.linspace(np.log(0.01), np.log(90), 100))
vec_maxc = np.exp(np.linspace(np.log(0.001), np.log(15), 100))

for i in range(len(vec_slope)-1):
    for j in range(len(vec_maxc)-1):
        stat = fn([0.5*(vec_slope[i]+vec_slope[i+1]), 0.5*(vec_maxc[j]+vec_maxc[j+1])])
        if np.isfinite(stat):
            stat_col_all = max(0.0001, min(0.9999, (stat - min(col_bounds)) / (max(col_bounds) - min(col_bounds))))
            col = cmap_cus(stat_col_all)
        else:
            col = 'tab:gray'

        ax.fill_between(x=[vec_slope[i], vec_slope[i+1]], y1=[vec_maxc[j], vec_maxc[j]],
                        y2=[vec_maxc[j+1], vec_maxc[j+1]], facecolor=col)

ax.set_ylim((0.001, 15))
ax.set_xlim((0.1, 90))
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel('Maximum absolute curvature (10$^{2}$ m$^{-1}$)')
ax.set_xlabel('Slope (degrees)')
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")

ax.text(-0.025, 1.025, 'b', transform=ax.transAxes, ha='right', va='bottom',fontweight='bold', fontsize=14)

# Save to file
plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_S8_final.png', dpi=400)