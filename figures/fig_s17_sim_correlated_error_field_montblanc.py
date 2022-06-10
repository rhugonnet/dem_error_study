"""Plotting of Figure S17: simulated correlated error fields around the Mont-Blanc summit"""
import gstools as gs
import matplotlib.pyplot as plt
import numpy as np
from geoutils import Raster
import xdem
import time
import pandas as pd
import matplotlib.colors as colors
import cartopy.crs as ccrs
import matplotlib.patches as mpatches

# Open data
fn_dem = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/Mont-Blanc_2017-10-25_DEM_5m.tif'
n_sim = 1
fn_hs = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/Mont-Blanc_2017-10-25_DEM_5m_hillshade.tif'
r = Raster(fn_dem)

# Around Mont-Blanc
crop_ext = [333500, 5076000, 335500, 5078000]
r.crop(crop_ext)
hs = Raster(fn_hs)
hs.crop(crop_ext)

fn_hetsce = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/df_heteroscedas_slope_curv.csv'
fn_vgm = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/df_vgm_std_sta.csv'

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
# Grid/Raster of 1000 x 1000 pixels
x = np.arange(0, shape[0]) * r.res[0]
y = np.arange(0, shape[1]) * r.res[0]

_, params = xdem.spatialstats.fit_sum_model_variogram(list_model=['Gau', 'Gau', 'Gau'], empirical_variogram=df_v,
                                                         bounds=[(0, 200), (0, 9), (500, 5000), (0, 9), (2000, 15000), (0,9)],
                                                         p0=[100, 1.5, 2000,1.5, 5000,1.5])

fn = xdem.spatialstats.interp_nd_binning(df_h, list_var_names=['slope_mid', 'maxc_mid'], statistic='nmad', min_count=30)
maxabsc[maxabsc>50] = 50
dh_err = fn((slope.data, maxabsc))

nmad_stable = 1.60
# nmad_stable = np.nanmedian(dh_err)

# Need to specify the rescale factor to match skgstat and gstools
model_s_alone = gs.Gaussian(dim=2, var=1, len_scale=params[0], rescale=2)

model_s = gs.Gaussian(dim=2, var=params[1], len_scale=params[0], rescale=2)
model_l = gs.Gaussian(dim=2, var=params[3], len_scale=params[2], rescale=2)
model_l2 = gs.Gaussian(dim=2, var=params[5], len_scale=params[4], rescale=2)

sim_slope_dems, sim_aspect_dems = (np.empty((6, n_sim,) + shape, dtype=np.float32) for i in range(2))

i=1
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
noisy_stationary_lr_dem = dem + nmad_stable * field_s + + nmad_stable * (field_l + field_l2)

# Heteroscedasticity with correlated noise
# !! Careful !! The long-range noise is scaled to the average variance, as it is not linked to heteroscedasticity
noisy_hetsce_sr_dem = dem + dh_err * field_s_alone
noisy_hetsce_lr_dem = dem + dh_err * field_s + nmad_stable * (field_l + field_l2)

# Function to plot a submap
def add_submap(fig, slices_grid, array, cmap, col_bounds, label, pos_colorbar=None, add_colorbar=True, label_colorbar=None, add_panel_letter=None, add_scale=False):

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
    # ax0.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
    ax0.text(0.5, 1.025, label, transform=ax0.transAxes, ha='center', va='bottom', fontweight='bold', fontsize=9, zorder=20)

    if add_panel_letter is not None:
        ax0.text(0.05, 0.95, add_panel_letter, transform=ax0.transAxes, ha='left', va='top', fontweight='bold', fontsize=14,
                zorder=20)

    if pos_colorbar is None:
        # pos = [0.2, -0.15, 0.6, 0.05]
        pos = [1.05, 0.2, 0.05, 0.6]
    else:
        pos = pos_colorbar

    if add_scale:
        y_extent = tmp_disp.bounds.top - tmp_disp.bounds.bottom
        x_extent = tmp_disp.bounds.right - tmp_disp.bounds.left
        ax0.add_patch(mpatches.Rectangle((crop_ext[2] - x_extent / 20 - 1000, crop_ext[1] + y_extent / 5), 500, 75,
                                        edgecolor='black', facecolor='black', transform=ccrs.UTM(32), zorder=10,
                                        linewidth=0.5))
        ax0.add_patch(mpatches.Rectangle((crop_ext[2] - x_extent / 20 - 500, crop_ext[1] + y_extent / 5), 500, 75,
                                        edgecolor='black', facecolor='white', transform=ccrs.UTM(32), zorder=10,
                                        linewidth=0.5))
        ax0.text(crop_ext[2] - x_extent / 20 - 1000, crop_ext[1] + y_extent / 5 - 20, '0', ha='center', va='top',
                transform=ccrs.UTM(32), zorder=10)
        ax0.text(crop_ext[2] - x_extent / 20 - 500, crop_ext[1] + y_extent / 5 - 20, '0.5', ha='center', va='top',
                transform=ccrs.UTM(32), zorder=10)
        ax0.text(crop_ext[2] - x_extent / 20 - 500, crop_ext[1] + y_extent / 5 - 150, 'km', ha='center', va='top',
                transform=ccrs.UTM(32), zorder=10)
        ax0.text(crop_ext[2] - x_extent / 20 - 0, crop_ext[1] + y_extent / 5 - 20, '1', ha='center', va='top',
                transform=ccrs.UTM(32), zorder=10)

    if add_colorbar:
        cbaxes = ax0.inset_axes(pos, zorder=10)

        norm = colors.Normalize(vmin=min(col_bounds), vmax=max(col_bounds))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cb = plt.colorbar(sm, cax=cbaxes, ticks=[-2, -1, 0, 1, 2], orientation='vertical', extend='both', shrink=0.2)
        cb.set_label(label_colorbar)


# Initiate figure and plot submaps
fig = plt.figure(figsize=(6, 12))
grid = plt.GridSpec(40, 20, wspace=0.1, hspace=0.1)

col_bounds = np.array([-2., -0.5, 0., 0.5,  2.])
cmap = plt.get_cmap('RdYlBu')
add_submap(fig, slices_grid=(slice(1, 9), slice(0, 10)), array=pixel_noise * nmad_stable, cmap=cmap, col_bounds=col_bounds,
           label='Homosc.,\nno corr.', add_colorbar=False, add_panel_letter='a', add_scale=True)
add_submap(fig, slices_grid=(slice(11, 19), slice(0, 10)), array=nmad_stable * field_s_alone, cmap=cmap, col_bounds=col_bounds,
           label='Homosc.,\nshort-range', add_colorbar=False, add_panel_letter='c')

add_submap(fig, slices_grid=(slice(31, 39), slice(0, 10)), array=nmad_stable * field_s + nmad_stable * field_l, cmap=cmap, col_bounds=col_bounds,
           label='Homosc.,\nlong-range', add_colorbar=False, add_panel_letter='g')

add_submap(fig, slices_grid=(slice(1, 9), slice(10, 20)), array=pixel_noise * dh_err, cmap=cmap, col_bounds=col_bounds,
           label='Heterosc.,\nno corr.', add_colorbar=False, add_panel_letter='b')
add_submap(fig, slices_grid=(slice(11, 19), slice(10, 20)), array=dh_err * field_s_alone, cmap=cmap, col_bounds=col_bounds,
           label='Heterosc.,\nshort-range', label_colorbar='Simulated random elevation error (m)', pos_colorbar=[1.05, -0.7, 0.075, 1.4], add_panel_letter='d')
add_submap(fig, slices_grid=(slice(31, 39), slice(10, 20)), array=dh_err * field_s + dh_err * field_l, cmap=cmap, col_bounds=col_bounds,
           label='Heterosc.,\nlong-range', add_colorbar=False, add_panel_letter='h')

add_submap(fig, slices_grid=(slice(21, 29), slice(0, 10)), array=nmad_stable * field_l, cmap=cmap, col_bounds=col_bounds,
           label='Homosc.,\nonly long-range', add_colorbar=False, add_panel_letter='e')
add_submap(fig, slices_grid=(slice(21, 29), slice(10, 20)), array=dh_err * field_l, cmap=cmap, col_bounds=col_bounds,
           label='Heterosc.,\nonly long-range', add_colorbar=False, add_panel_letter='f')

# Save to file
plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_S17_final.png', dpi=300)