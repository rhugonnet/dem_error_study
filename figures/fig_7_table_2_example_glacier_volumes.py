"""Plotting of Figure 7: example of propagation of errors to glacier volume changes for the Mont-Blanc case study"""
import os
from typing import Callable
from scipy.spatial.distance import pdist
import numpy as np
import pandas as pd
from geoutils import Raster, Vector
import geoutils as gu
import xdem
import matplotlib.pyplot as plt
import seaborn as sns
import skgstat
import pyproj

fn_dem = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/Mont-Blanc_2017-10-25_DEM_5m.tif'
fn_dh = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/dh_NK_Deramp_final.tif'
r = Raster(fn_dem)
dh_r = Raster(fn_dh)
dh = gu.spatial_tools.get_array_and_mask(dh_r)[0]
fn_shp_glacier = '/home/atom/data/inventory_products/RGI/00_rgi60_neighb_merged/11_rgi60_CentralEurope/11_rgi60_CentralEurope.shp'
glacier_inventory = Vector(fn_shp_glacier)
glacier_inventory.crop2raster(r)

# Open estimates of heteroscedasticity and spatial correlation of errors
fn_hetsce = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/df_heteroscedas_slope_curv.csv'
fn_vgm = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/df_vgm_std_sta.csv'

# Model those
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


_, params_lr = xdem.spatialstats.fit_sum_model_variogram(list_model=['Gau', 'Sph', 'Sph'], empirical_variogram=df_v,
                                                         bounds=[(0, 200), (0, 9), (500, 5000), (0, 9), (2000, 15000), (0,9)],
                                                         p0=[100, 1.5, 2000,1.5, 5000,1.5])

# _, params_llr = xdem.spatialstats.fit_sum_model_variogram(list_model=['Sph', 'Sph', 'Sph'], empirical_variogram=df_v[df_v.bins<28000],
#                                                       bounds=[(0, 50), (0, 1), (100, 5000), (0, 1), (1000, 25000), (0, 1)],
#                                                       p0=[25, 0.5, 2500, 0.5, 10000, 0.5])

fn = xdem.spatialstats.interp_nd_binning(df_h, list_var_names=['slope_mid', 'maxc_mid'], statistic='nmad', min_count=30)
maxabsc[maxabsc>50] = 50
dh_err = fn((slope.data, maxabsc))
dh_err[dh_err>np.nanpercentile(dh_err,99)] = np.nanpercentile(dh_err,99)
nmad_stable = 1.60

list_name = ['Homosc., no corr.', 'Homosc., short-range', 'Homosc., long-range',
             'Heterosc., no corr.', 'Heterosc., short-range', 'Heterosc., long-range', 'Heterosc., supra-long']

list_df = []
list_rgiids = glacier_inventory.ds['RGIId'].values
list_df_propag = []
for rgiid in list_rgiids:

    gla_ds = glacier_inventory.ds[glacier_inventory.ds['RGIId'].values==rgiid]
    gla_shp = Vector(gla_ds)
    lat, lon = gla_ds['CenLat'].values[0], gla_ds['CenLon'].values[0]
    gla_mask = gla_shp.rasterize(r, in_value=1) == 1

    # Let's only use glaciers with at least 70% data
    glacier_pixels = np.count_nonzero(gla_mask)
    if glacier_pixels == 0:
        continue
    valid_glacier_pixels = np.count_nonzero(np.isfinite(dh[gla_mask]))
    perc_valid_data = valid_glacier_pixels/glacier_pixels * 100.

    mean_slope = np.nanmean(slope[gla_mask])

    dh_gla = np.nanmean(dh[gla_mask])

    if perc_valid_data > 85.:

        # Glacier area
        area = glacier_pixels*r.res[0]*r.res[1]

        # Error with stationary variance and no correlation
        err_dh_stationary_nocorr = nmad_stable / np.sqrt(glacier_pixels)

        # Error with non-stationary variance and no correlation
        err_dh_hetsce_nocorr = np.nanmean(dh_err[gla_mask]) / np.sqrt(glacier_pixels)

        # Error with stationary variance and short-range correlation
        neff_sr = xdem.spatialstats.neff_circ(area=area, list_vgm=[(params_lr[0]/2, 'Gau', params_lr[1])])
        err_dh_stationary_corr_sr = nmad_stable / np.sqrt(neff_sr)

        # Error with non-stationary variance and short-range correlation
        err_dh_hetsce_corr_sr = np.nanmean(dh_err[gla_mask]) / np.sqrt(neff_sr)

        # Error with stationary variance and long-range correlations
        # !! Careful !! Need to scale by 1/2 for the Gaussian model to match that of scikit-gstat
        neff_lr = xdem.spatialstats.neff_circ(area=area, list_vgm=[(params_lr[0]/2, 'Gau', params_lr[1]),
                                                                   (params_lr[2], 'Sph', params_lr[3]),
                                                                   (params_lr[4], 'Sph', params_lr[5])])

        neff_supralr = xdem.spatialstats.neff_circ(area=area, list_vgm=[(params_lr[0]/2, 'Gau', params_lr[1]),
                                                                   (params_lr[2], 'Sph', params_lr[3]),
                                                                   (20000, 'Sph', 1.5*params_lr[5])])

        neff_lr_only = xdem.spatialstats.neff_circ(area=area, list_vgm=[(params_lr[2], 'Sph', params_lr[3])])
        neff_lr2_only = xdem.spatialstats.neff_circ(area=area, list_vgm=[(params_lr[4], 'Sph', params_lr[5])])

        err_dh_stationary_corr_lr = nmad_stable / np.sqrt(neff_lr)

        # Error with non-stationary variance and short-range correlation
        err_dh_hetsce_corr_lr = np.nanmean(dh_err[gla_mask]) / np.sqrt(neff_lr)

        err_dh_hetsce_corr_supralr = np.nanmean(dh_err[gla_mask]) / np.sqrt(neff_supralr)

        list_err = [err_dh_stationary_nocorr, err_dh_stationary_corr_sr, err_dh_stationary_corr_lr,
                    err_dh_hetsce_nocorr, err_dh_hetsce_corr_sr, err_dh_hetsce_corr_lr, err_dh_hetsce_corr_supralr]

        list_df_cat = []
        for i in range(len(list_err)):
            df_cat = pd.DataFrame()
            df_cat = df_cat.assign(rgiid = [rgiid], area = [area], mean_slope = [mean_slope], dh=[dh_gla], err_dh=[list_err[i]], err_category=[list_name[i]])
            list_df_cat.append(df_cat)

        df_sub = pd.concat(list_df_cat)

        df_propag = pd.DataFrame()
        df_propag = df_propag.assign(rgiid = [rgiid], area = [area], lat=[lat], lon=[lon],
                                     std_err_sr=[params_lr[1]/np.sqrt(neff_sr)], std_err_lr=[params_lr[3]/np.sqrt(neff_lr_only)],
                                     std_err_lr2=[params_lr[5]/np.sqrt(neff_lr2_only)])

        list_df.append(df_sub)
        list_df_propag.append(df_propag)

df = pd.concat(list_df)
df_p = pd.concat(list_df_propag)

# Text values on underestimation for a glacier of 10 km²
area_example = 10*10**6

# Random: divide by the square root of the pixel number
err_random_gla = nmad_stable / np.sqrt(area_example/(5**2))
neff_example_sr = xdem.spatialstats.neff_circ(area=area_example, list_vgm=[(params_lr[0] / 2, 'Gau', params_lr[1])])

err_short_range = nmad_stable / np.sqrt(neff_example_sr)

neff_example_lr = xdem.spatialstats.neff_circ(area=area_example, list_vgm=[(params_lr[0]/2, 'Gau', params_lr[1]),
                                                                   (params_lr[2], 'Sph', params_lr[3]),
                                                                   (params_lr[4], 'Sph', params_lr[5])])
err_long_range = nmad_stable / np.sqrt(neff_example_lr)

random_underestimation_factor = err_long_range/err_random_gla
short_range_underestimation_factor = err_long_range/err_short_range

print('Glacier of '+str(area_example)+' km², error underestimated by:')
print(str(short_range_underestimation_factor)+ ' using short range vario')
print(str(random_underestimation_factor)+' using no corr.')

def double_sum_covar(coords: np.ndarray, areas: np.ndarray, errors: list[np.ndarray],
                     vgm_funcs: list[Callable]):
    """
    Double sum of covariance for euclidean coordinates
    :param coords: Spatial support (typically, pixel) coordinates
    :param areas:  Area of supports
    :param errors: Standard errors of supports
    :param vgm_funcs: Variogram function
    :return:
    """

    n = len(coords)
    pds = pdist(coords)
    var = 0
    for i in range(n):
        for j in range(n):

            # For index calculation of the pairwise distance, see https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html
            if i == j:
                d = 0
            elif i < j:
                ind = n * i + j - ((i + 2) * (i + 1)) // 2
                d = pds[ind]
            else:
                ind = n * j + i - ((j + 2) * (j + 1)) // 2
                d = pds[ind]

            for k in range(len(vgm_funcs)):
                var += errors[k][i] * errors[k][j] * (1 - vgm_funcs[k](d)) * areas[i] * areas[j]

    total_area = sum(areas)
    se_dsc = np.sqrt(var / total_area ** 2)

    return se_dsc


def vgm_func_short(h):

    return skgstat.models.gaussian(h, params_lr[0], params_lr[1]) / params_lr[1]


def vgm_func_long(h):

    return skgstat.models.spherical(h, params_lr[2], params_lr[3]) / params_lr[3]

def vgm_func_long2(h):
    return skgstat.models.spherical(h, params_lr[4], params_lr[5]) / params_lr[5]


xx, yy = gu.projtools.reproject_from_latlon((df_p.lat.values, df_p.lon.values), out_crs=pyproj.CRS(32632))

massif_err_nocorr = np.sqrt(np.sum(df_p.std_err_sr**2 * df_p.area**2) / np.sum(df_p.area)**2)
coords = np.array((xx, yy)).T
massif_err_shortrange = double_sum_covar(coords=coords, areas=df_p.area.values, errors=[df_p.std_err_sr.values], vgm_funcs=[vgm_func_short])
massif_err_longrange = double_sum_covar(coords=coords, areas=df_p.area.values, errors=[df_p.std_err_sr.values, df_p.std_err_lr.values, df_p.std_err_lr2.values],
                                        vgm_funcs=[vgm_func_short, vgm_func_long, vgm_func_long2])
area_massif = np.sum(df_p.area) / 1000000

rgiid_mdg_bossons = ['RGI60-11.03646', 'RGI60-11.03647']

df_mdg_bossons = df_p[df_p.rgiid.isin(rgiid_mdg_bossons)]
mdg_bossons_err_nocorr = np.sqrt(np.sum(df_mdg_bossons.std_err_sr**2 * df_mdg_bossons.area**2) / np.sum(df_mdg_bossons.area)**2)
xx, yy = gu.projtools.reproject_from_latlon((df_mdg_bossons.lat.values, df_mdg_bossons.lon.values), out_crs=pyproj.CRS(32632))
coords = np.array((xx, yy)).T
mdg_bossons_err_shortrange = double_sum_covar(coords=coords, areas=df_mdg_bossons.area.values, errors=[df_mdg_bossons.std_err_sr.values], vgm_funcs=[vgm_func_short])
mdg_bossons_err_longrange = double_sum_covar(coords=coords, areas=df_mdg_bossons.area.values, errors=[df_mdg_bossons.std_err_sr.values, df_mdg_bossons.std_err_lr.values, df_mdg_bossons.std_err_lr2.values],
                                        vgm_funcs=[vgm_func_short, vgm_func_long, vgm_func_long2])
area_bossons_taconnaz = np.sum(df_mdg_bossons.area) / 1000000

rgiid_griaz_bourgeat = ['RGI60-11.03280', 'RGI60-11.03290']


df_griaz_bourgeat = df_p[df_p.rgiid.isin(rgiid_griaz_bourgeat)]
griaz_bourgeat_err_nocorr = np.sqrt(np.sum(df_griaz_bourgeat.std_err_sr**2 * df_griaz_bourgeat.area**2) / np.sum(df_griaz_bourgeat.area)**2)
xx, yy = gu.projtools.reproject_from_latlon((df_griaz_bourgeat.lat.values, df_griaz_bourgeat.lon.values), out_crs=pyproj.CRS(32632))
coords = np.array((xx, yy)).T
griaz_bourgeat_err_shortrange = double_sum_covar(coords=coords, areas=df_griaz_bourgeat.area.values, errors=[df_griaz_bourgeat.std_err_sr.values], vgm_funcs=[vgm_func_short])
griaz_bourgeat_err_longrange = double_sum_covar(coords=coords, areas=df_griaz_bourgeat.area.values, errors=[df_griaz_bourgeat.std_err_sr.values, df_griaz_bourgeat.std_err_lr.values, df_griaz_bourgeat.std_err_lr2.values],
                                        vgm_funcs=[vgm_func_short, vgm_func_long, vgm_func_long2])
area_griaz_bourgeat = np.sum(df_griaz_bourgeat.area) / 1000000

# Print table:
print('Values (No corr., Short range, Long range)')
print('Massif:')
print([massif_err_nocorr, massif_err_shortrange, massif_err_longrange])
print('Bosson Taconnaz:')
print([mdg_bossons_err_nocorr, mdg_bossons_err_shortrange, mdg_bossons_err_longrange])
print('Griaz Bourgeat:')
print([griaz_bourgeat_err_nocorr, griaz_bourgeat_err_shortrange, griaz_bourgeat_err_longrange])

# Initiate figure

fig = plt.figure(figsize=(6.5, 10.5))
grid = plt.GridSpec(22, 22, wspace=0.1, hspace=0.1)

# First, an horizontal axis on top to plot the sample histograms


bins_area = [0.01, 0.1, 1, 5, 30]
list_gla_nb = []
df['Glacier area category (km²)'] = ''
for i in range(len(bins_area)-1):
    ind = np.logical_and(df.area/1000000 >= bins_area[i], df.area/1000000 < bins_area[i+1])
    df.loc[ind, ['Glacier area category (km²)']] = str(bins_area[i])+'–'+str(bins_area[i+1])
    nb_gla = len(np.unique(df[ind].rgiid.values))
    list_gla_nb.append(nb_gla)
df = df.sort_values(by='Glacier area category (km²)')

ax = fig.add_subplot(grid[0:2, :10])

for i, gla_nb in enumerate(list_gla_nb):
    ax.fill_between([i-1/5, i+1/5], [0] * 2, [gla_nb] * 2,
                     facecolor='black', alpha=1,
                     edgecolor='black', linewidth=0.5)

ax.vlines(x=np.arange(0.5, len(list_gla_nb)), ymin=0, ymax=max(list_gla_nb)*1.1, linestyles='dotted', colors='grey', linewidth=0.5)
ax.set_xlim((-0.5, len(list_gla_nb)-0.5))
ax.set_ylim((0, max(list_gla_nb)*1.1))
ax.set_xticklabels([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylabel('Number of\nglaciers')

ax = fig.add_subplot(grid[2:10, :10])

sns.boxplot(ax=ax, x="Glacier area category (km²)", y="err_dh", hue="err_category", hue_order=list_name[0:6],
                 data=df, palette={list_name[0]:'white', list_name[1]:'lightgrey', list_name[2]:'darkgrey' ,
                                   list_name[3]:'lightgreen', list_name[4]:'green', list_name[5]:'darkgreen'},
                 fliersize=0, linewidth=1)

ax.vlines(x=np.arange(0.5, len(list_gla_nb)), ymin=-1, ymax=5, linestyles='dotted', colors='grey', linewidth=0.5)
ax.legend(ncol=1)
ax.set_ylabel('Uncertainty in the\nmean elevation change (1$\sigma$, m)')
ax.set_yscale('log')
ax.set_ylim((0.00001, 2))
ax.set_xlim((-0.5, len(list_gla_nb)-0.5))
ax.text(0.03, 0.95, 'a', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14,
         zorder=20)

# Same with glacier mean slope

bins_slope = [15, 25, 35, 45, 55]
list_gla_nb = []
df['Glacier mean slope category (km²)'] = ''
for i in range(len(bins_area)-1):
    ind = np.logical_and(df.mean_slope >= bins_slope[i], df.mean_slope < bins_slope[i+1])
    df.loc[ind, ['Glacier mean slope category (km²)']] = str(bins_slope[i])+'-'+str(bins_slope[i+1])
    nb_gla = len(np.unique(df[ind].rgiid.values))
    list_gla_nb.append(nb_gla)
df = df.sort_values(by='Glacier mean slope category (km²)')


ax = fig.add_subplot(grid[12:14, :10])

for i, gla_nb in enumerate(list_gla_nb):
    ax.fill_between([i-1/5, i+1/5], [0] * 2, [gla_nb] * 2,
                     facecolor='black', alpha=1,
                     edgecolor='black', linewidth=0.5)

ax.vlines(x=np.arange(0.5, len(list_gla_nb)), ymin=0, ymax=max(list_gla_nb)*1.1, linestyles='dotted', colors='grey', linewidth=0.5)
ax.set_xlim((-0.5, len(list_gla_nb)-0.5))
ax.set_ylim((0, max(list_gla_nb)*1.1))
ax.set_xticklabels([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylabel('Number of\nglaciers')

ax = fig.add_subplot(grid[14:, :10])

sns.boxplot(ax=ax, x="Glacier mean slope category (km²)", y="err_dh", hue="err_category", hue_order=list_name[0:6],
                 data=df, palette={list_name[0]:'white', list_name[1]:'lightgrey', list_name[2]:'darkgrey' ,
                                   list_name[3]:'lightgreen', list_name[4]:'green', list_name[5]:'darkgreen'},
                 fliersize=0, linewidth=1)

ax.vlines(x=np.arange(0.5, len(list_gla_nb)), ymin=-1, ymax=5, linestyles='dotted', colors='grey', linewidth=0.5)
l =ax.legend()
l.remove()
ax.set_ylabel('Uncertainty in the\nmean elevation change (1$\sigma$, m)')
ax.set_yscale('log')
ax.set_ylim((0.001, 2))
ax.set_xlim((-0.5, len(list_gla_nb)-0.5))
ax.text(0.03, 0.95, 'b', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14,
         zorder=20)

# LAST PANEL: per-glacier uncertainty intersections to zero


list_valid_rgiids = np.unique(df.rgiid)

df_s_nc = df[df.err_category=='Homosc., no corr.']
df_s_nc = df_s_nc.sort_values(by='area')
df_s_nc['glacier_nb'] = np.arange(0, len(list_valid_rgiids))

df_h_csr = df[df.err_category=='Heterosc., short-range']
df_h_csr = df_h_csr.sort_values(by='area')
df_h_csr['glacier_nb'] = np.arange(0, len(list_valid_rgiids))

df_h_clr = df[df.err_category=='Heterosc., long-range']
df_h_clr = df_h_clr.sort_values(by='area')
df_h_clr['glacier_nb'] = np.arange(0, len(list_valid_rgiids))

df_h_supralr = df[df.err_category=='Heterosc., supra-long']
df_h_supralr = df_h_supralr.sort_values(by='area')
df_h_supralr['glacier_nb'] = np.arange(0, len(list_valid_rgiids))

ax = fig.add_subplot(grid[:, -3:])

areas = df_h_csr.area.values

for i in range(len(list_valid_rgiids)):

    ax.fill_between([0, areas[i]/1000000], [i+-0.3] * 2, [i+0.3] * 2, facecolor='black',
        edgecolor='black', linewidth=0.5)

ax.set_yticks([])
ax.set_ylim((-0.5, len(list_valid_rgiids)))
ax.set_xlabel('Area of\nglacier (km²)')
ax.set_xscale('log')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax = fig.add_subplot(grid[:, 13:-3])

ax.vlines(x=0, ymin=0, ymax=len(list_valid_rgiids), colors='black', linestyles='dashed')

for i in range(len(list_valid_rgiids)):
    ax.fill_between(
        [df_h_clr.dh.values[i] - 2 * df_h_clr.err_dh.values[i], df_h_clr.dh.values[i] + 2 * df_h_clr.err_dh.values[i]],
        [df_h_clr.glacier_nb.values[i] - 0.4] * 2, [df_h_clr.glacier_nb.values[i] + 0.4] * 2, facecolor='darkgreen',
        edgecolor='black', linewidth=0.5, alpha=0.8)
    ax.fill_between(
        [df_h_csr.dh.values[i] - 2 * df_h_csr.err_dh.values[i], df_h_csr.dh.values[i] + 2 * df_h_csr.err_dh.values[i]],
        [df_h_csr.glacier_nb.values[i] - 0.4] * 2, [df_h_csr.glacier_nb.values[i] + 0.4] * 2, facecolor='lightgrey',
        edgecolor='black', linewidth=0.5, alpha=0.8)
    ax.fill_between(
        [df_s_nc.dh.values[i] - 5 * df_s_nc.err_dh.values[i], df_s_nc.dh.values[i] + 5 * df_s_nc.err_dh.values[i]],
        [df_s_nc.glacier_nb.values[i] - 0.4] * 2, [df_s_nc.glacier_nb.values[i] + 0.4] * 2, facecolor='white',
        edgecolor='black', linewidth=0.5, alpha=0.8)

inters_zero_s_nc = np.count_nonzero(np.logical_and(df_s_nc.dh.values - 2 * df_s_nc.err_dh.values < 0, df_s_nc.dh.values + 2 * df_s_nc.err_dh.values > 0))
inters_zero_h_csr = np.count_nonzero(np.logical_and(df_h_csr.dh.values - 2 * df_h_csr.err_dh.values < 0, df_h_csr.dh.values + 2 * df_h_csr.err_dh.values > 0))
inters_zero_h_clr = np.count_nonzero(np.logical_and(df_h_clr.dh.values - 2 * df_h_clr.err_dh.values < 0, df_h_clr.dh.values + 2 * df_h_clr.err_dh.values > 0))

df_h_supralr_large = df_h_supralr[df_h_supralr.area>200000.]
inters_zero_h_supralr_large = np.count_nonzero(np.logical_and(df_h_supralr_large.dh.values - 2 * df_h_supralr_large.err_dh.values < 0,
                                                              df_h_supralr_large.dh.values + 2 * df_h_supralr_large.err_dh.values > 0))

df_h_supralr_small = df_h_supralr[df_h_supralr.area<200000.]
inters_zero_h_supralr_small = np.count_nonzero(np.logical_and(df_h_supralr_small.dh.values - 2 * df_h_supralr_small.err_dh.values < 0,
                                                              df_h_supralr_small.dh.values + 2 * df_h_supralr_small.err_dh.values > 0))

perc_s_nc = inters_zero_s_nc / len(list_valid_rgiids) * 100.
perc_h_csr = inters_zero_h_csr / len(list_valid_rgiids) * 100.
perc_h_clr = inters_zero_h_clr / len(list_valid_rgiids) * 100.

perc_h_supralr_large = inters_zero_h_supralr_large / len(df_h_supralr_large) * 100.
perc_h_supralr_small = inters_zero_h_supralr_small / len(df_h_supralr_small) * 100.

print('Percentage of 95% ranges intersecting zeros for:')
print('No correlation: {:.2f}'.format(perc_s_nc))
print('Short-range only: {:.2f}'.format(perc_h_csr))
print('Short- and long-range: {:.2f}'.format(perc_h_clr))
print('Short- and long-range, and swath range for big glaciers (>0.2 km²): {:.2f}'.format(perc_h_supralr_large))
print('Short- and long-range, and swath range for small glaciers (<0.2 km²): {:.2f}'.format(perc_h_supralr_small))

ax.scatter(x=df_s_nc.dh.values, y=df_s_nc.glacier_nb, marker='o', facecolor='black', s=5)

ax.set_xlim((-6, 6))
ax.set_ylabel('Individual glaciers ordered by increasing area')

ax.set_xlabel('Mean elevation\nchange (m)')

p = ax.fill_between([],[],[], facecolor='white', edgecolor='black', linewidth=0.5,
                    label='Uncertainty\n(2$\sigma$) for\nhomosc. and\nno corr.: \n{:.0f}% intersect\nzero'.format(perc_s_nc))
p2 = ax.fill_between([],[],[], facecolor='lightgrey', edgecolor='black', linewidth=0.5,
                     label='Uncertainty\n(2$\sigma$) for\nheterosc. and\nshort-range:\n{:.0f}% intersect\nzero'.format(perc_h_csr))
p3= ax.fill_between([],[],[], facecolor='darkgreen', edgecolor='black', linewidth=0.5,
                    label='Uncertainty\n(2$\sigma$) for\nheterosc. and\nlong-range:\n{:.0f}% intersect\nzero'.format(perc_h_clr))
p4 = ax.scatter([], [], facecolor='black', marker='o', s=5)
ax.legend(loc='lower right', bbox_to_anchor=(1.85, 0.45), framealpha=0.9)
ax.set_ylim((-0.5, len(list_valid_rgiids)))

ax.set_yticks(np.arange(0, len(list_valid_rgiids)))
labels_glaciers = ['' for i in range(len(list_valid_rgiids))]
labels_glaciers[0] = '1$^{st}$'
labels_glaciers[-1] = str(len(list_valid_rgiids))+'$^{th}$'
for i in np.arange(10, 90, 10):
    labels_glaciers[i-1] = str(i)+'$^{th}$'
ax.set_yticklabels(labels_glaciers)
ax.text(0.05, 0.98, 'c', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14,
         zorder=20)

plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_7_final.png', dpi=400)
