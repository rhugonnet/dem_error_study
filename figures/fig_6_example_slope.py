"""Plotting of Figure 6: example of propagation of errors to terrain slope and aspect for the Mont-Blanc case study"""
import os
import gstools as gs
import matplotlib.pyplot as plt
import numpy as np
from geoutils import Raster
import xdem
import time
import pandas as pd
import seaborn as sns
import matplotlib.colors as colors
import cartopy.crs as ccrs

# Open file
fn_dem = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/Mont-Blanc_2017-10-25_DEM_5m.tif'
n_sim = 200
fn_hs = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/Mont-Blanc_2017-10-25_DEM_5m_hillshade.tif'
r = Raster(fn_dem)

# Crop around Mont-Blanc
crop_ext = [333500, 5076000, 335500, 5078000]
r.crop(crop_ext)
hs = Raster(fn_hs)
hs.crop(crop_ext)

# Open files with estimates of heteroscedasticity and spatial correlation
fn_hetsce = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/df_heteroscedas_slope_curv.csv'
fn_vgm = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/df_vgm_std_sta.csv'

# Model those estimates
df_h = pd.read_csv(fn_hetsce)
df_v = pd.read_csv(fn_vgm)
df_v = df_v[df_v.bins<30000]
df_v.err_exp /= np.sqrt(100)/2
std_fac = np.nanmean(df_v.exp.values[-3:])
df_v.exp /= std_fac
df_v.err_exp /= std_fac
dem = np.copy(r.data.data).squeeze()
dem[dem==r.nodata] = np.nan

slope, aspect, planc, profc = xdem.terrain.get_terrain_attribute(dem, resolution=r.res[0], attribute=['slope', 'aspect', 'planform_curvature',
                                                                                'profile_curvature'])
maxabsc = np.maximum(np.abs(planc), np.abs(profc))

shape = np.shape(dem)

# Grid definition
x = np.arange(0, shape[0]) * r.res[0]
y = np.arange(0, shape[1]) * r.res[0]

# Fit sum of Gaussian models (only a tiny difference with Spherical for long range, works better with Fourier simulation)
fun1, params = xdem.spatialstats.fit_sum_model_variogram(list_model=['Gau', 'Gau', 'Gau'], empirical_variogram=df_v,
                                                         bounds=[(0, 200), (0, 9), (500, 5000), (0, 9), (2000, 15000), (0,9)],
                                                         p0=[100, 1.5, 2000,1.5, 5000,1.5])

fn = xdem.spatialstats.interp_nd_binning(df_h, list_var_names=['slope_mid', 'maxc_mid'], statistic='nmad', min_count=30)
maxabsc[maxabsc>50] = 50
dh_err = fn((slope.data, maxabsc))

# NMAD on complete stable terrain (not only cropped)
nmad_stable = 1.60
# nmad_stable = np.nanmedian(dh_err)

# Need to specify the rescale factor to match Gaussian parameters between skgstat and gstools
model_s_alone = gs.Gaussian(dim=2, var=1, len_scale=params[0], rescale=2)

model_s = gs.Gaussian(dim=2, var=params[1], len_scale=params[0], rescale=2)
model_l = gs.Gaussian(dim=2, var=params[3], len_scale=params[2], rescale=2)
model_l2 = gs.Gaussian(dim=2, var=params[5], len_scale=params[4], rescale=2)

sim_slope_dems, sim_aspect_dems = (np.empty((6, n_sim,) + shape, dtype=np.float32) for i in range(2))

for i in range(n_sim):

    print('Working on simulation '+str(i+1))

    print('Generating random field...')

    t0 = time.time()

    # Using GSTools, let's generate a correlated signal at two different length: 5 and 100 (spherical)
    srf_s_alone = gs.SRF(model_s_alone, mode_no=100)
    srf_s = gs.SRF(model_s, mode_no=100)
    srf_l = gs.SRF(model_l, mode_no=100)
    srf_l2 = gs.SRF(model_l2, mode_no=100)

    # We combine the two random correlated fields (e.g, short-range could represent resolution, and long-range the noise)
    field_s_alone = srf_s_alone.structured([x, y])

    field_s = srf_s((x, y), mesh_type='structured')
    field_l = srf_l((x, y), mesh_type='structured')
    field_l2 = srf_l2((x, y), mesh_type='structured')

    # Stationary variance with purely random noise
    pixel_noise = np.random.normal(0, 1, size=np.shape(dem))
    noisy_stationary_dem = dem + pixel_noise * nmad_stable

    # Heteroscedasticity with purely random noise
    noisy_hetsce_dem = dem + pixel_noise * dh_err

    # Stationary variance with correlated noise (short, and short+long range)
    noisy_stationary_sr_dem = dem + nmad_stable * field_s_alone
    noisy_stationary_lr_dem = dem + nmad_stable * field_s + nmad_stable * (field_l + field_l2)

    # Heteroscedasticity with correlated noise
    # !! Careful !! The long-range noise is scaled to the average variance, as it is not linked to heteroscedasticity
    noisy_hetsce_sr_dem = dem + dh_err * field_s_alone
    noisy_hetsce_lr_dem = dem + dh_err * field_s + nmad_stable * (field_l + field_l2)

    t1 = time.time()

    print('Elapsed: {:.1f} seconds'.format(t1-t0))

    print('Deriving slopes...')

    # Derive attribute for each simulation
    slope_stationary, aspect_stationary = xdem.terrain.get_terrain_attribute(noisy_stationary_dem, resolution=r.res[0], attribute=['slope', 'aspect'])
    slope_hetsce, aspect_hetsce = xdem.terrain.get_terrain_attribute(noisy_hetsce_dem, resolution=r.res[0], attribute=['slope', 'aspect'])
    slope_stationary_sr, aspect_stationary_sr = xdem.terrain.get_terrain_attribute(noisy_stationary_sr_dem, resolution=r.res[0], attribute=['slope', 'aspect'])
    slope_stationary_lr, aspect_stationary_lr = xdem.terrain.get_terrain_attribute(noisy_stationary_lr_dem, resolution=r.res[0], attribute=['slope', 'aspect'])
    slope_hetsce_sr, aspect_hetsce_sr = xdem.terrain.get_terrain_attribute(noisy_hetsce_sr_dem, resolution=r.res[0], attribute=['slope', 'aspect'])
    slope_hetsce_lr, aspect_hetsce_lr = xdem.terrain.get_terrain_attribute(noisy_hetsce_lr_dem, resolution=r.res[0], attribute=['slope', 'aspect'])

    t2 = time.time()
    print('Elapsed: {:.1f} seconds'.format(t2-t1))

    # Save in 4D array
    sim_slope_dems[0, i, :, :] = slope_stationary
    sim_slope_dems[1, i, :, :] = slope_hetsce
    sim_slope_dems[2, i, :, :] = slope_stationary_sr
    sim_slope_dems[3, i, :, :] = slope_stationary_lr
    sim_slope_dems[4, i, :, :] = slope_hetsce_sr
    sim_slope_dems[5, i, :, :] = slope_hetsce_lr

    sim_aspect_dems[0, i, :, :] = aspect_stationary
    sim_aspect_dems[1, i, :, :] = aspect_hetsce
    sim_aspect_dems[2, i, :, :] = aspect_stationary_sr
    sim_aspect_dems[3, i, :, :] = aspect_stationary_lr
    sim_aspect_dems[4, i, :, :] = aspect_hetsce_sr
    sim_aspect_dems[5, i, :, :] = aspect_hetsce_lr


# Define bins to compute statistics on
bins_slope = [0, 5, 10, 15, 20, 30, 40, 50, 70, 90]
bins_curv = [0, 0.2, 0.5, 1, 2, 3, 4, 6, 10, 20, 50]

northness = np.cos(sim_aspect_dems * np.pi / 180)
eastness = np.sin(sim_aspect_dems * np.pi / 180)

list_slope_map, list_maxnortheast_map, list_df_bp, list_df_bp_northeast = ([] for i in range(4))
run_names = ['stationary_random', 'hetsce_random', 'stationary_shortrange', 'stationary_longrange', 'hetsce_shortrange', 'hetsce_longrange']
for i in range(6):


    slope_1sig = (np.nanpercentile(sim_slope_dems[i, :, :, :], 84, axis=0)\
                        - np.nanpercentile(sim_slope_dems[i, :, :, :], 16, axis=0)) / 2
    northness_1sig = (np.nanpercentile(northness[i, :, :, :], 84, axis=0)\
                        - np.nanpercentile(northness[i, :, :, :], 16, axis=0)) / 2
    eastness_1sig = (np.nanpercentile(eastness[i, :, :, :], 84, axis=0)\
                        - np.nanpercentile(eastness[i, :, :, :], 16, axis=0)) / 2
    maxnortheast_1sig = np.maximum(northness_1sig, eastness_1sig)

    for j in range(len(bins_slope) - 1):
        # Subset by slope category
        subset = np.logical_and(slope >= bins_slope[j], slope < bins_slope[j + 1])
        sub_slope = slope_1sig[subset]
        sub_northeast = maxnortheast_1sig[subset]
        # Store in dataframe
        df_subset = pd.DataFrame()
        df_subset = df_subset.assign(err_slope=sub_slope, run=[run_names[i]]*len(sub_slope))
        bin_name = str(bins_slope[j]) + '–' + str(bins_slope[j + 1])
        df_subset['slope_category'] = bin_name
        list_df_bp.append(df_subset)

        df_subset_northeast = pd.DataFrame()
        df_subset_northeast = df_subset_northeast.assign(err_northeast=sub_northeast, run=[run_names[i]] * len(sub_slope))
        bin_name = str(bins_slope[j]) + '–' + str(bins_slope[j + 1])
        df_subset_northeast['slope_category'] = bin_name
        list_df_bp_northeast.append(df_subset_northeast)

    list_slope_map.append(slope_1sig)
    list_maxnortheast_map.append(maxnortheast_1sig)

list_bins_slope_2 = []
for j in range(len(bins_slope) - 1):
    list_bins_slope_2.append(str(bins_slope[j]) + '–\n' + str(bins_slope[j + 1]))

df_bp = pd.concat(list_df_bp)
df_bp_northeast = pd.concat(list_df_bp_northeast)

# Subset for boxplot of uncertainty with slope categories
orig_names = ['stationary_random', 'stationary_longrange','hetsce_random', 'hetsce_longrange']
df_bp_sub = df_bp[df_bp.run.isin(orig_names)]
df_bp_northeast_sub = df_bp_northeast[df_bp_northeast.run.isin(orig_names)]

names = ['Homosc., no corr.', 'Homosc., long-range', 'Heterosc., no corr.', 'Heterosc., long-range']
for i, oname in enumerate(orig_names):
    df_bp_sub.loc[df_bp_sub.run == oname, 'run'] = names[i]
    df_bp_northeast_sub.loc[df_bp_northeast_sub.run == oname, 'run'] = names[i]


# Define a submap function

def add_submap(fig, slices_grid, array, cmap, col_bounds, label, pos_colorbar=None, add_colorbar=True, label_colorbar=None, add_panel_letter=None):

    ax0 = fig.add_subplot(grid[slices_grid[0], slices_grid[1]], projection=ccrs.UTM(32), label=label)

    tmp_disp = r.copy()
    ext = [tmp_disp.bounds[0], tmp_disp.bounds[2], tmp_disp.bounds[1], tmp_disp.bounds[3]]

    tmp_disp.data[0, :, :] = array

    cb = []
    cb_val = np.linspace(0, 1, len(col_bounds))
    for j in range(len(cb_val)):
        cb.append(cmap(cb_val[j]))
    cmap_cus2 = colors.LinearSegmentedColormap.from_list('my_cb', list(
        zip((col_bounds - min(col_bounds)) / (max(col_bounds - min(col_bounds))), cb)), N=1000)
    cmap_cus2.set_bad(color='None')
    ax0.imshow(tmp_disp.data[0, :, :], extent=ext, transform=ccrs.UTM(32), vmin=min(col_bounds), vmax=max(col_bounds), cmap=cmap_cus2,
              interpolation=None, zorder=3, alpha=0.85)
    ax0.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
    ax0.text(0.5, 1.025, label, transform=ax0.transAxes, ha='center', va='bottom', fontweight='bold', fontsize=9, zorder=20)

    if add_panel_letter is not None:
        ax0.text(-0.1, 1.2, add_panel_letter, transform=ax0.transAxes, ha='left', va='top', fontweight='bold', fontsize=14,
                zorder=20)

    if pos_colorbar is None:
        pos = [-0.5, -0.15, 1, 0.05]
        # pos = [1.05, 0.2, 0.05, 0.6]
    else:
        pos = pos_colorbar

    if add_colorbar:
        cbaxes = ax0.inset_axes(pos, zorder=10)

        norm = colors.Normalize(vmin=min(col_bounds), vmax=max(col_bounds))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cb = plt.colorbar(sm, cax=cbaxes, ticks=col_bounds, orientation='horizontal', extend='both', shrink=0.2)
        cb.set_label(label_colorbar)


# Start figure
fig = plt.figure(figsize=(7, 11))
grid = plt.GridSpec(40, 35, wspace=0.1, hspace=0.1)

# First, an horizontal axis on top to plot the sample histograms

ax = fig.add_subplot(grid[:6, 18:])

list_nb_pixel = []
for i in range(len(bins_slope)-1):
    ind_pixel = np.logical_and(df_bp.run.values == 'stationary_random', df_bp.slope_category== str(bins_slope[i])+'–'+str(bins_slope[i+1]))
    nb_pixel = np.count_nonzero(ind_pixel)
    list_nb_pixel.append(nb_pixel)
    ax.fill_between([i-0.3, i+0.3], [0]*2, [nb_pixel], facecolor='black')

ax.vlines(np.arange(0.5, len(bins_slope)-1), ymin=-5, ymax=np.max(list_nb_pixel)*1.1, colors='tab:gray', linestyles='dashed', linewidths=0.75)

ax.set_xticks([])
ax.set_ylabel('Sample count')
ax.set_ylim((100, np.max(list_nb_pixel)*1.1))
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_yscale('log')
ax.set_xlim((-0.5, len(bins_slope)-1.5))
ax.text(0.025, 0.9, 'c', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14,
         zorder=20)

# Boxplots of slope error depending on slope
ax = fig.add_subplot(grid[6:23, 18:])

sns.boxplot(ax=ax, x="slope_category", y="err_slope", hue="run", hue_order=names,
                 data=df_bp_sub, palette={names[0]:'white', names[1]:'darkgrey', names[2]:'lightgreen' ,names[3]:'darkgreen'},
                 fliersize=0, linewidth=1)
ax.vlines(np.arange(0.5, len(bins_slope)-1), ymin=-5, ymax=40, colors='tab:gray', linestyles='dashed', linewidths=0.75)

ax.set_ylim((-0.5, 17.25))
ax.set_xlabel('Slope categories (degrees)')
ax.set_ylabel('Uncertainty in slope (1$\sigma$, degrees)')
ax.legend(loc='upper right')
ax.set_xticklabels([])
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
ax.set_xlim((-0.5, len(bins_slope)-1.5))

# Boxplots of aspect error depending on slope
ax = fig.add_subplot(grid[23:, 18:])

sns.boxplot(ax=ax, x="slope_category", y="err_northeast", hue="run", hue_order=names,
                 data=df_bp_northeast_sub, palette={names[0]:'white', names[1]:'darkgrey', names[2]:'lightgreen' ,names[3]:'darkgreen'},
                 fliersize=0, linewidth=1)
ax.vlines(np.arange(0.5, len(bins_slope)-1), ymin=-1, ymax=2, colors='tab:gray', linestyles='dashed', linewidths=0.5)
ax.set_ylim((-0.05, 1.05))
l = ax.legend()
l.remove()

# ax.text(0.025, 0.96, 'd', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14,
#          zorder=20)
ax.set_xlabel('Slope categories (degrees)')
ax.set_ylabel('Maximum of uncertainty in northness or eastness (1$\sigma$)')
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
ax.set_xticklabels(list_bins_slope_2)
ax.set_xlim((-0.5, len(bins_slope)-1.5))

# Finally, the minimaps
color1 = colors.to_rgba('black')
color2 = colors.to_rgba('white')
cmap = colors.LinearSegmentedColormap.from_list('my_cmap_hs', [color1, color2], 256)
col_bounds = np.array([0., 127.5, 255.])
add_submap(fig, slices_grid=(slice(1, 9), slice(0, 8)), array=hs.data[0, :, :], cmap=cmap, col_bounds=col_bounds, label='Hillshade',
           add_colorbar=False, add_panel_letter='a')

cmap = plt.get_cmap('Reds', 100)
col_bounds = np.array([0., 20., 40., 60., 80.])
add_submap(fig, slices_grid=(slice(1, 9), slice(8, 16)), array=slope, cmap=cmap, col_bounds=col_bounds, label='Slope',
           label_colorbar='Slope (degrees)')

col_bounds = np.array([-1., 0.,  1.])
cmap = plt.get_cmap('RdBu')
add_submap(fig, slices_grid=(slice(11, 19), slice(0, 8)), array=np.cos(aspect * np.pi / 180), cmap=cmap, col_bounds=col_bounds,
           label='Northness', add_colorbar=False)
add_submap(fig, slices_grid=(slice(11, 19), slice(8, 16)), array=np.sin(aspect * np.pi / 180), cmap=cmap, col_bounds=col_bounds,
           label_colorbar='Northess or Eastness', label='Eastness')

col_bounds = np.array([0., 5, 10.])
cmap = plt.get_cmap('Purples')
add_submap(fig, slices_grid=(slice(22, 30), slice(0, 8)), array=list_slope_map[0], cmap=cmap, col_bounds=col_bounds,
           label='Homosc.,\nno corr.', add_colorbar=False, add_panel_letter='b')

add_submap(fig, slices_grid=(slice(22, 30), slice(8, 16)), array=list_slope_map[5], cmap=cmap, col_bounds=col_bounds,
           label_colorbar='Uncertainty in slope (1$\sigma$, degrees)', label='Heterosc.,\nlong-range')

col_bounds = np.array([0., 0.25,  0.5])
cmap = plt.get_cmap('Greens')
add_submap(fig, slices_grid=(slice(33, 41), slice(0, 8)), array=list_maxnortheast_map[0], cmap=cmap, col_bounds=col_bounds,
           label='Homosc.,\nno corr.', add_colorbar=False)

add_submap(fig, slices_grid=(slice(33, 41), slice(8, 16)), array=list_maxnortheast_map[5], cmap=cmap, col_bounds=col_bounds,
           label_colorbar='Maximum uncertainty in\nnorthness or eastness (1$\sigma$)', label='Heterosc.,\nlong-range')

# Save to file
plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_6_final.png', dpi=400)
