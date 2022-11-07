"""Plotting of Figure S18: assymetry of slope and aspect errors for the Mont-Blanc summit"""
import gstools as gs
import matplotlib.pyplot as plt
import numpy as np
from geoutils import Raster
import xdem
import time
import pandas as pd

# Open data
fn_dem = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/Pleiades_Mont-Blanc_2017-10-25_DEM_5m.tif'
n_sim = 200
r = Raster(fn_dem)

# Crop around Mont-Blanc
crop_ext = [334500, 5077000, 335000, 5077500]
r.crop(crop_ext)

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

    slope_stationary, aspect_stationary = xdem.terrain.get_terrain_attribute(noisy_stationary_dem, resolution=r.res[0], attribute=['slope', 'aspect'])
    slope_hetsce, aspect_hetsce = xdem.terrain.get_terrain_attribute(noisy_hetsce_dem, resolution=r.res[0], attribute=['slope', 'aspect'])
    slope_stationary_sr, aspect_stationary_sr = xdem.terrain.get_terrain_attribute(noisy_stationary_sr_dem, resolution=r.res[0], attribute=['slope', 'aspect'])
    slope_stationary_lr, aspect_stationary_lr = xdem.terrain.get_terrain_attribute(noisy_stationary_lr_dem, resolution=r.res[0], attribute=['slope', 'aspect'])
    slope_hetsce_sr, aspect_hetsce_sr = xdem.terrain.get_terrain_attribute(noisy_hetsce_sr_dem, resolution=r.res[0], attribute=['slope', 'aspect'])
    slope_hetsce_lr, aspect_hetsce_lr = xdem.terrain.get_terrain_attribute(noisy_hetsce_lr_dem, resolution=r.res[0], attribute=['slope', 'aspect'])

    t2 = time.time()
    print('Elapsed: {:.1f} seconds'.format(t2-t1))

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



northness = np.cos(sim_aspect_dems * np.pi / 180)
eastness = np.sin(sim_aspect_dems * np.pi / 180)

## Initiate figure
fig = plt.figure(figsize=(6.5, 6.5))

grid = plt.GridSpec(22, 23, wspace=0.1, hspace=0.1)

ax = fig.add_subplot(grid[0:10, 0:10])

# Slope error distribution on low slopes
ind = np.logical_and(slope>0, slope<20)

statio = sim_slope_dems[0, :, ind].flatten()
hetsce_lr = sim_slope_dems[5, :, ind].flatten()
ax.hist(statio, alpha=0.5, facecolor='white', bins=50, edgecolor='tab:gray', linewidth=0.5, density=True, label='Homosc.,\nno corr.')
ax.hist(hetsce_lr, alpha=0.5, facecolor='tab:green', bins=50, edgecolor='tab:gray', linewidth=0.5, density=True, label='Heterosc.,\nlong-range')
ax.text(0.5, 0.8, '0° < Initial slope < 20°', transform=ax.transAxes, ha='center', va='top')
ax.text(0.05, 0.95, 'a', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14,
         zorder=20)
ax.set_ylabel('Density')
ax.set_xlabel('Simulated slope (degrees)')
ax.legend(loc='center right')

# Slope error distribution on high slopes
ax = fig.add_subplot(grid[0:10, 13:])

ind = np.logical_and(slope>70, slope<90)

statio = sim_slope_dems[0, :, ind].flatten()
hetsce_lr = sim_slope_dems[5, :, ind].flatten()
ax.hist(statio, alpha=0.5, facecolor='white', bins=50, edgecolor='tab:gray', linewidth=0.5, density=True)
ax.hist(hetsce_lr, alpha=0.5, facecolor='tab:green', bins=50, edgecolor='tab:gray', linewidth=0.5, density=True)
ax.text(0.5, 0.8, '70° < Initial slope < 90°', transform=ax.transAxes, ha='center', va='top')
ax.text(0.05, 0.95, 'b', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14,
         zorder=20)
ax.set_ylabel('Density')
ax.set_ylim((0, 0.16))
ax.set_xlabel('Simulated slope (degrees)')

# Aspect error distribution on low slopes
ax = fig.add_subplot(grid[12:, 0:10])

ind = np.logical_and(slope>0, slope<20)

statio = northness[0, :, ind].flatten()
hetsce_lr = northness[5, :, ind].flatten()
ax.hist(statio, alpha=0.5, facecolor='white', bins=50, edgecolor='tab:gray', linewidth=0.5, density=True)
ax.hist(hetsce_lr, alpha=0.5, facecolor='tab:green', bins=50, edgecolor='tab:gray', linewidth=0.5, density=True)
ax.text(0.5, 0.8, '0° < Initial slope < 20°', transform=ax.transAxes, ha='center', va='top')
ax.text(0.05, 0.95, 'c', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14,
         zorder=20)
ax.set_ylabel('Density')
ax.set_xlabel('Simulated northness')

# Aspect error distribution on high slopes
ax = fig.add_subplot(grid[12:, 13:])

ind = np.logical_and(slope>70, slope<90)

statio = northness[0, :, ind].flatten()
hetsce_lr = northness[5, :, ind].flatten()
ax.hist(statio, alpha=0.5, facecolor='white', bins=50, edgecolor='tab:gray', linewidth=0.5, density=True)
ax.hist(hetsce_lr, alpha=0.5, facecolor='tab:green', bins=50, edgecolor='tab:gray', linewidth=0.5, density=True)
ax.text(0.5, 0.8, '70° < Initial slope < 90°', transform=ax.transAxes, ha='center', va='top')
ax.text(0.05, 0.95, 'd', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14,
         zorder=20)
ax.set_ylabel('Density')
ax.set_xlabel('Simulated northness')

# Save to file
plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_S18_final.png', dpi=400)