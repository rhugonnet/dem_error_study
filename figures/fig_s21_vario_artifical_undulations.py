"""Plotting of Figure S21: artificial undulations to constrain correlated errors in swath direction"""
import xdem
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import skgstat as skg
from matplotlib.patches import Patch
from matplotlib.legend_handler import HandlerPatch

# Simulate undulations scenario 1
ampli = 1
res = 1
freq = 10

nx, ny = (100, 100)
shape = (nx, ny)
x = np.arange(nx)
y = np.arange(ny)
xx, yy = np.meshgrid(x, y)

along_angle = 0

altrack = -xx * np.sin(np.deg2rad(along_angle)) + yy * np.cos(np.deg2rad(along_angle))
dh_simu = ampli * np.sin(altrack * res * 2 * np.pi / freq)

coords = np.dstack((xx.flatten(), yy.flatten())).squeeze()
vals = dh_simu.flatten()

# Estimate variogram
V = skg.Variogram(coords, vals, n_lags=100)
bins, exps = V.get_empirical()
counts = V.bin_count

df = pd.DataFrame()
df = df.assign(bins=bins, exp=exps, counts=counts, err_exp=np.nan*np.ones(len(bins)))
df = df[df.bins<=115]

fun, params = xdem.spatialstats.fit_sum_model_variogram(['Sph'], df)
dfr = df.copy()
dfr.exp = df.exp.rolling(5, min_periods=1).min().drop_duplicates()
dfr[np.logical_and(dfr.bins>1,dfr.bins<=5)] = np.nan
fun_double, params_double = xdem.spatialstats.fit_sum_model_variogram(['Sph', 'Sph'], dfr)

# Simulate undulations scenario 2
freq = 30

altrack2 = -xx * np.sin(np.deg2rad(along_angle)) + yy * np.cos(np.deg2rad(along_angle))
dh_simu2 = ampli * np.sin(altrack2 * res * 2 * np.pi / freq)
vals2 = dh_simu2.flatten()

# Estimate variogram
V2 = skg.Variogram(coords, vals2, n_lags=100)
bins2, exps2 = V2.get_empirical()
counts2 = V2.bin_count

df2 = pd.DataFrame()
df2 = df2.assign(bins=bins2, exp=exps2, counts=counts2, err_exp=np.nan*np.ones(len(bins)))
df2 = df2[df2.bins<=115]

fun2, params2 = xdem.spatialstats.fit_sum_model_variogram(['Sph'], df2)

df2r = df2.copy()
df2r.exp = df2.exp.rolling(15, min_periods=1).min().drop_duplicates()
df2r[np.logical_and(df2r.bins>1,df2r.bins<=30)] = np.nan
fun2_double, params2_double = xdem.spatialstats.fit_sum_model_variogram(['Sph', 'Sph'], df2r)

xmin = 0
xmax = max(df.bins)

# Initiate figure
fig = plt.figure(figsize=(12, 12))
grid = plt.GridSpec(22, 22, wspace=0.1, hspace=0.1)

# Plot artificial undulations, first scenario
ax = fig.add_subplot(grid[0:10, 0:10])

ax.imshow(dh_simu, cmap='RdYlBu', vmin=-1, vmax=1)
ax.text(-0.025, 1.05, 'a', transform=ax.transAxes, ha='right', va='bottom', fontweight='bold', fontsize=14)

# Plot artificial undulations, second scenario
ax = fig.add_subplot(grid[12:, 0:10])
ax.text(-0.025, 1.05, 'b', transform=ax.transAxes, ha='right', va='bottom', fontweight='bold', fontsize=14)

ax.imshow(dh_simu2, cmap='RdYlBu')

# Plot variograms
ax0 = fig.add_subplot(grid[:3, 12:])

ax0.set_xticks([])

# Plot the histogram manually with fill_between
interval_var = [0] + list(df.bins)
for i in range(len(df)):
    width = interval_var[i + 1] - interval_var[i]
    mid = interval_var[i] + width / 2
    count = df['counts'].values[i]
    ax0.fill_between([mid - width / 2, mid + width /2], [0] * 2, [count] * 2,
                     facecolor='black', alpha=1,
                     edgecolor='none', linewidth=0.5)
    # ax0.vlines(mid - width / 2, ymin=[0], ymax=1.2 * max(df['counts'].values), colors='tab:gray',
    #            linestyles='dashed', linewidths=0.5)
ax0.set_ylabel('Pairwise\nsample\ncount')
# Scientific format to avoid undesired additional space on the label side
ax0.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

# Ignore warnings for log scales
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((0, 1.2 * max(df['counts'].values)))
ax0.text(0.05, 0.95, 'c', transform=ax0.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)

ax1 = fig.add_subplot(grid[3:, 12:])

bins_model = np.linspace(0, 150, 1000)

bins_center = np.subtract(df.bins, np.diff([0] + df.bins.tolist()) / 2)

bins_center_r = np.subtract(dfr.bins, np.diff([0] + dfr.bins.tolist()) / 2)
bins_center_r2 = np.subtract(df2r.bins, np.diff([0] + df2r.bins.tolist()) / 2)

ax1.scatter(bins_center, df.exp, marker='x', color='tab:blue')
ax1.scatter(bins_center_r, dfr.exp, marker='x', s=80, color='tab:blue', linewidths=3)
ax1.scatter(bins_center, df2.exp, marker='x', color='tab:orange')
ax1.scatter(bins_center_r2, df2r.exp, marker='x', s=80, color='tab:orange', linewidths=3)
ax1.plot(bins_model, fun(bins_model), color='tab:blue')
ax1.plot(bins_model, fun_double(bins_model), color='tab:blue', linestyle='dashed')
ax1.plot(bins_model, fun2(bins_model), color='tab:orange')
ax1.plot(bins_model, fun2_double(bins_model), color='tab:orange', linestyle='dashed')
ax1.hlines(0.5, xmin=xmin, xmax=xmax, linestyles='dotted', colors='black', label='Global variance')
ax1.set_xlim((xmin, xmax))
ax1.set_ylabel('Variance')
ax1.set_xlabel('Spatial lag')

ax1.scatter([], [], marker='x', label='Empirical variogram', color='black')
ax1.scatter([], [], marker='x', s=80, linewidths=3, label='Rolling minimum values', color='black')
ax1.plot([], [], linestyle='solid', label='Double range fit on all values', color='black')
ax1.plot([], [], linestyle='dashed', label='Double range fit on minimum values', color='black')

handles, labels = ax1.get_legend_handles_labels()
order = [3, 4, 0, 1, 2]
l0 = ax1.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='lower right')
# l0.set_zorder(30)

p0 = Patch(facecolor='tab:blue', edgecolor='None', label='Variogram of panel a')
p1 = Patch(facecolor='tab:orange', edgecolor='None', label='Variogram of panel b')

ax1.legend([p0, p1], ['Variogram of panel a', 'Variogram of panel b'],
           handler_map={p0 : HandlerPatch(), p1: HandlerPatch()},
           framealpha=1, loc='upper right')
# l1.set_zorder(30)
ax1.add_artist(l0)

# ax1.legend(loc='lower right')

# Save to file
plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_S21_final.png', dpi=300)