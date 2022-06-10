"""Plotting of Figure S13: improve pairwise sampling for variogram estimation on grid data"""
import numpy as np
import skgstat as skg
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Define random state
rnd = np.random.RandomState(np.random.MT19937(np.random.SeedSequence(42)))

shape = (500, 500)
x = np.arange(0, shape[0])
y = np.arange(0, shape[1])
xx, yy = np.meshgrid(x, y)
vals = np.random.normal(0, 1, size=shape)

# Flatten everything because we don't care about the 2D at this point
coords = np.dstack((xx.flatten(), yy.flatten())).squeeze()

# Completely random subsetting
rnd = np.random.RandomState(np.random.MT19937(np.random.SeedSequence(42)))
subset = rnd.choice(len(coords), 200, replace=False)
coords_sub = coords[subset, :]

extent=(x[0], x[-1], y[0], y[-1])
ratio_subsample = 0.2
samples = int(200 / 3.72)
res = np.mean([(extent[1] - extent[0])/(shape[0]-1),(extent[3] - extent[2])/(shape[1]-1)])
V_r = skg.Variogram(coordinates=coords_sub, values=vals.flatten()[subset], normalize=False)

# Disk equidistant subsetting
rnd = np.random.RandomState(np.random.MT19937(np.random.SeedSequence(42)))
rems_mp = skg.RasterEquidistantMetricSpace(coords, shape=shape, extent=extent, samples=samples, rnd=rnd, runs=1, verbose=True)
V = skg.Variogram(rems_mp, values=vals.flatten(), normalize=False)
coords_eq = coords[rems_mp.eqidx[0]]
coords_center = coords[rems_mp.cidx[0]]
rnd = np.random.RandomState(np.random.MT19937(np.random.SeedSequence(42)))
idx_center = rnd.choice(len(coords), size=1, replace=False)
center = rems_mp._centers[0]
radius = np.sqrt(1. / ratio_subsample * samples / np.pi) * res

equidistant_radii = [0.]
increasing_rad = radius
max_dist = np.sqrt((extent[1] - extent[0])**2 + (extent[3] - extent[2])**2)
while increasing_rad < max_dist:
    equidistant_radii.append(increasing_rad)
    increasing_rad *= np.sqrt(2)
equidistant_radii.append(max_dist)

V_r2 = skg.Variogram(coordinates=coords_sub, values=vals.flatten()[subset], normalize=False, bin_func=equidistant_radii[1:])
rnd = np.random.RandomState(np.random.MT19937(np.random.SeedSequence(42)))
rems_mp = skg.RasterEquidistantMetricSpace(coords, shape=shape, extent=extent, samples=samples, rnd=rnd, runs=1, verbose=True)
V2 = skg.Variogram(rems_mp, values=vals.flatten(), normalize=False, bin_func=equidistant_radii[1:])
V2_bins = [np.sum(np.logical_and(V2.distance>=equidistant_radii[i], V2.distance<equidistant_radii[i+1])) for i in range(len(equidistant_radii)-1)]

# Plot figure
fig = plt.figure(figsize=(12, 12))

grid = plt.GridSpec(22, 24, wspace=0.1, hspace=0.1)

# First, just a grid of random coordinates
ax = fig.add_subplot(grid[0:10, 0:11])

ax.scatter(coords_sub[:, 0], coords_sub[:, 1], marker='o', color='tab:grey', s=10, label='Full grid subset')
ax.set_xlim((0, shape[0]))
ax.set_ylim((0, shape[1]))
ax.text(0.025, 0.975, 'a', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14,
         zorder=20)

ax.legend(loc='lower center')

# Then, a grid of coordinates sampled for pairwise comparison between disk and rings
ax = fig.add_subplot(grid[0:10, 13:])


ax.scatter(coords_eq[:, 0], coords_eq[:, 1], marker='o', color='tab:grey', s=10, label='Ring subsets')
ax.scatter(coords_center[:, 0], coords_center[:, 1], marker='o', color='tab:orange', s=10, label='Disk subset')

ax.add_patch(mpatches.Circle(center, radius=equidistant_radii[1], edgecolor='black', facecolor='None', linewidth=1, linestyle='dashed'))
for r in equidistant_radii[2:]:
    ax.add_patch(mpatches.Circle(center, radius=r, edgecolor='black', facecolor='None', linewidth=1, linestyle='dashed'))

ax.scatter(center[0], center[1], marker='x', s=50, label='Center', color='black')
ax.text(0.025, 0.975, 'b', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14,
         zorder=20)

ax.legend(loc='lower center')

ax.set_xlim((0, shape[0]))
ax.set_ylim((0, shape[1]))

# Pairwise count with spatial lag for the random sampling
ax = fig.add_subplot(grid[12:16, 0:11])

bin_edges = [0] + list(V_r.bins)
bin_count = V_r.bin_count
for i in range(len(bin_edges)-1):
    ax.fill_between([bin_edges[i], bin_edges[i+1]], y1=[0]*2, y2=[bin_count[i]]*2, edgecolor='white', facecolor='tab:grey')
ax.set_ylim((0, 1.25*max(bin_count)))
ax.set_ylabel('Pairwise\nsample count')
ax.set_xlabel('Distance (linear)')
ax.text(0.5, 0.95, 'Total pairs: {:,}'.format(np.sum(V_r.bin_count)), transform=ax.transAxes
        , ha='center', va='top')
ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
ax.text(0.025, 0.95, 'c', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14,
         zorder=20)

# Same with log X-axis
ax = fig.add_subplot(grid[18:22, 0:11])

bin_edges = [7] + list(V_r2.bins)
bin_count = V_r2.bin_count
for i in range(len(bin_edges)-1):
    ax.fill_between([bin_edges[i], bin_edges[i+1]], y1=[0]*2, y2=[bin_count[i]]*2, edgecolor='white', facecolor='tab:grey')
ax.set_ylim((0, 1.25*max(bin_count)))
ax.set_ylabel('Pairwise\nsample count')
ax.set_xlabel('Distance (logarithmic)')
ax.set_xscale('log')
ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
ax.text(0.025, 0.95, 'e', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14,
         zorder=20)

# Pairwise count with spatial lag for the disk/ring sampling
ax = fig.add_subplot(grid[12:16, 13:])

bin_edges = [0] + list(V.bins)
bin_count = V.bin_count
for i in range(len(bin_edges)-1):
    ax.fill_between([bin_edges[i], bin_edges[i+1]], y1=[0]*2, y2=[bin_count[i]]*2, edgecolor='white', facecolor='tab:grey')
ax.set_ylim((0, 1.25*max(bin_count)))
ax.set_ylabel('Pairwise\nsample count')
ax.set_xlabel('Distance (linear)')
ax.text(0.5, 0.95, 'Total pairs: {:,}'.format(np.sum(V.bin_count)), transform=ax.transAxes
        ,ha='center', va='top')
ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
ax.text(0.025, 0.95, 'd', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14,
         zorder=20)

# Same with log X-axis
ax = fig.add_subplot(grid[18:22, 13:])

bin_edges = [7] + list(V2.bins)
bin_count = V2.bin_count
for i in range(len(bin_edges)-1):
    ax.fill_between([bin_edges[i], bin_edges[i+1]], y1=[0]*2, y2=[bin_count[i]]*2, edgecolor='white', facecolor='tab:grey')
ax.set_ylim((0, 1.25*max(bin_count)))
ax.set_ylabel('Pairwise\nsample count')
ax.set_xlabel('Distance (logarithmic)')
ax.set_xscale('log')
ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
ax.text(0.025, 0.95, 'f', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14,
         zorder=20)

# Save to file
plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_S13_final.png', dpi=400)