"""Plotting of Figure S19: impact of short-range correlation close to a 3x3 kernel size on slope errors"""
import gstools as gs
import matplotlib.pyplot as plt
import numpy as np
from geoutils import Raster
import xdem
import time
import pandas as pd
import seaborn as sns

# Open data
fn_dem = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/Pleiades_Mont-Blanc_2017-10-25_DEM_5m.tif'
n_sim = 200
r = Raster(fn_dem)

# Crop around Mont-Blanc
crop_ext = [333500, 5076000, 335500, 5078000]
r.crop(crop_ext)

dem = np.copy(r.data.data).squeeze()
dem[dem==r.nodata] = np.nan

slope, aspect, planc, profc = xdem.terrain.get_terrain_attribute(dem, resolution=r.res[0], attribute=['slope', 'aspect', 'planform_curvature',
                                                                                'profile_curvature'])
maxabsc = np.maximum(np.abs(planc), np.abs(profc))

shape = np.shape(dem)
# Grid/Raster of 1000 x 1000 pixels
x = np.arange(0, shape[0]) * r.res[0]
y = np.arange(0, shape[1]) * r.res[0]


nmad_stable = 1.60
# nmad_stable = np.nanmedian(dh_err)

# Need to specify the rescale factor to match skgstat and gstools

list_len_scale=[1.25, 2.5, 5, 10, 15, 25, 50, 100]
list_model_s = []
for len_scale in list_len_scale:
    model_s = gs.Gaussian(dim=2, var=1, len_scale=len_scale, rescale=2)
    list_model_s.append(model_s)

sim_slope_dems, sim_aspect_dems = (np.empty((8, n_sim,) + shape, dtype=np.float32) for i in range(2))

for i in range(n_sim):

    print('Working on simulation '+str(i+1))

    print('Generating random field...')

    t0 = time.time()

    for j in range(len(list_len_scale)):
        # Using GSTools, let's generate a correlated signal at two different length: 5 and 100 (spherical)
        srf_s_alone = gs.SRF(list_model_s[j], mode_no=100)

        # We combine the two random correlated fields (e.g, short-range could represent resolution, and long-range the noise)
        field_s_alone = srf_s_alone.structured([x, y])

        # Stationary variance with correlated noise (short, and short+long range)
        noisy_stationary_sr_dem = dem + nmad_stable * field_s_alone

        t1 = time.time()

        print('Elapsed: {:.1f} seconds'.format(t1-t0))

        print('Deriving slopes...')

        slope_stationary_sr, aspect_stationary_sr = xdem.terrain.get_terrain_attribute(noisy_stationary_sr_dem, resolution=r.res[0], attribute=['slope', 'aspect'])

        t2 = time.time()
        print('Elapsed: {:.1f} seconds'.format(t2-t1))

        sim_slope_dems[j, i, :, :] = slope_stationary_sr

run_names = []
for i in range(len(list_len_scale)):
    if i<1:
        run_names.append('Short range: {:.2f}'.format(list_len_scale[i]/5.)+' pixels')
    elif i<3:
        run_names.append('Short range: {:.1f}'.format(list_len_scale[i]/5.)+' pixels')
    else:
        run_names.append('Short range: {:.0f}'.format(list_len_scale[i] / 5.) + ' pixels')

list_bin_edges=[(0, 20), (70, 90)]
list_df_bp = []

for i in range(len(list_len_scale)):


    slope_1sig = (np.nanpercentile(sim_slope_dems[i, :, :, :], 84, axis=0)\
                        - np.nanpercentile(sim_slope_dems[i, :, :, :], 16, axis=0)) / 2

    for j in range(2):
        # Subset by slope category
        subset = np.logical_and(slope >= list_bin_edges[j][0], slope < list_bin_edges[j][1])
        sub_slope = slope_1sig[subset]
        # Store in dataframe
        df_subset = pd.DataFrame()
        df_subset = df_subset.assign(err_slope=sub_slope, run=[run_names[i]]*len(sub_slope))
        bin_name = str(list_bin_edges[j][0]) + '–' + str(list_bin_edges[j][1])
        df_subset['slope_category'] = bin_name
        list_df_bp.append(df_subset)

df_bp = pd.concat(list_df_bp)


# Initiate figure
fig = plt.figure(figsize=(6, 6))

grid = plt.GridSpec(22, 23, wspace=0.1, hspace=0.1)

ax = fig.add_subplot(grid[:, :])

sns.boxplot(ax=ax, x="slope_category", y="err_slope", hue="run", hue_order=run_names,
                 data=df_bp, #palette={run_names[0]:'white', run_names[1]:'lightgrey', run_names[2]:'darkgrey' ,run_names[3]:'tab:green'},
                 fliersize=0, linewidth=1)
ax.vlines(0.5, ymin=-0.5, ymax=12, colors='tab:gray', linestyles='dashed', linewidths=0.75)
ax.set_xlim((-0.5, 1.5))
ax.set_ylim((-0.5, 12))
ax.set_xlabel('Slope categories (degrees)')
ax.set_ylabel('Uncertainty in slope (1$\sigma$, degrees)')
ax.legend(loc='upper right')

# Save to file
plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_S19_final.png', dpi=400)

