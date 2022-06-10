"""Estimate heteroscedasticity and spatial correlation of errors for the Mont-Blanc case study"""
import numpy as np
import pandas as pd
import xdem
import geoutils as gu

# Open files
fn_ddem = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/dh_NK_Deramp_final.tif'
fn_pleiades = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/Mont-Blanc_2017-10-25_DEM_5m.tif'
fn_shp = '/home/atom/data/inventory_products/RGI/00_rgi60_neighb_merged/11_rgi60_CentralEurope/11_rgi60_CentralEurope.shp'
fn_forest = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/outlines/ESA_CCI_forest_simplified_delainey.shp'

# Custom mask to see more terrain on Fig. 5 panel b
# fn_forest = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/outlines/forest_Mt-Blanc.shp'

# Files on server
# fn_ddem = '/data/icesat/travail_en_cours/romain/dem_error_study/final_run/dh_NK_Deramp_final.tif'
# fn_pleiades = '/data/icesat/travail_en_cours/romain/dem_error_study/final_run/Mont-Blanc_2017-10-25_DEM_5m.tif'
# fn_shp = '/data/icesat/travail_en_cours/romain/data/outlines/rgi60/00_rgi60_neighb_merged/11_rgi60_CentralEurope/11_rgi60_CentralEurope.shp'
# fn_forest = '/data/icesat/travail_en_cours/romain/dem_error_study/final_run/outlines/ESA_CCI_forest_simplified_delainey.shp'


pleia_ddem = gu.Raster(fn_ddem)
ref_dem = gu.Raster(fn_pleiades)
glaciers_outlines = gu.Vector(fn_shp)
forest_outlines = gu.Vector(fn_forest)
mask_glacier = glaciers_outlines.create_mask(pleia_ddem)
mask_forest = forest_outlines.create_mask(pleia_ddem)

# Remove forest, very large outliers
pleia_ddem.data[mask_forest] = np.nan
pleia_ddem.data[np.abs(pleia_ddem.data)>500] = np.nan

slope, planc, profc = xdem.terrain.get_terrain_attribute(ref_dem, attribute=['slope', 'planform_curvature',
                                                                                'profile_curvature'])
maxabsc = np.maximum(np.abs(planc), np.abs(profc))

# 0/ Filter large outliers per category of slope and curvature
bins_slope = [0, 2.5, 5, 10, 15, 20, 30, 40, 50, 70, 90]
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
    # Subset by curvature category
    subset = np.logical_and(maxabsc >= bins_curv[i], maxabsc < bins_curv[i + 1])
    dh_sub = pleia_ddem.data[subset]
    # Remove very large outliers of the category
    med_sub = np.nanmedian(dh_sub)
    nmad_sub = xdem.spatialstats.nmad(dh_sub)
    # Remove outliers
    pleia_ddem.data[np.logical_and(subset, np.abs(pleia_ddem.data-med_sub) > 7 * nmad_sub)] = np.nan


# Subsample on stable terrain
pleia_ddem_sta = pleia_ddem.data[~mask_glacier]
slope_sta = slope.data[~mask_glacier]
maxabsc_sta = maxabsc[~mask_glacier]

# 1/ Estimate heteroscedasticity with slope and maximum curvature
df_sub = xdem.spatialstats.nd_binning(pleia_ddem_sta, list_var=[slope_sta, maxabsc_sta], list_var_names=['slope', 'maxc'], list_var_bins=(bins_slope, bins_curv))
df_sub['slope_mid'] = pd.IntervalIndex(df_sub.slope).mid.values
df_sub['maxc_mid'] = pd.IntervalIndex(df_sub.maxc).mid.values

# Save binned estimates to file
# df_sub.to_csv('/data/icesat/travail_en_cours/romain/dem_error_study/final_run/df_heteroscedas_slope_curv.csv')

# Interpolate with N-D binning
fn = xdem.spatialstats.interp_nd_binning(df_sub, list_var_names=['slope', 'maxc'], statistic='nmad', min_count=30)

# Create an error map, filter outliers, and standardize elevation differences
maxabsc[maxabsc>50] = 50
dh_err = fn((slope.data, maxabsc))
std_dh = pleia_ddem.data.data/dh_err
std_dh[np.abs(std_dh)>7*xdem.spatialstats.nmad(std_dh)] = np.nan
std_dh /= xdem.spatialstats.nmad(std_dh)
std_r = pleia_ddem.copy(new_array=std_dh)

del dh_err, maxabsc, slope, planc, profc, slope_sta, maxabsc_sta, pleia_ddem_sta
# Save standardized elevation difference map to file
std_r.save('/data/icesat/travail_en_cours/romain/dem_error_study/final_run/dh_NK_Deramp_std.tif')
# This one is saved with a custom mask to see more terrain Fig. 5
# std_r.save('/home/atom/ongoing/work_stderr_dem/case_study_montblanc/dh_NK_Deramp_std_customforest.tif')

# Subsample on glacier terrain
pleia_ddem_gla = pleia_ddem.copy()
pleia_ddem_gla.data[~mask_glacier] = np.nan
pleia_ddem.data.data[mask_glacier] = np.nan

std_dh_gla = np.copy(std_dh)
std_dh_gla[~mask_glacier] = np.nan
std_dh[mask_glacier] = np.nan

# 2/ Run patches method with a fixed random state for validation
# areas_emp = [100 * 2 ** i for i in range(18)]
# list_stderr_empirical, list_stderr_empirical_gla = ([] for i in range(2))
# for area_emp in areas_emp:
#
#     print('Working on patches for area size: '+str(area_emp))
#
#     #  First, sample intensively circular patches of a given area, and derive the mean elevation differences
#     df_patches = xdem.spatialstats.patches_method(std_dh.squeeze(), gsd=pleia_ddem.res[0], area=area_emp, n_patches=1000, random_state=42,
#                                                   perc_min_valid=80.)
#     df_patches_gla = xdem.spatialstats.patches_method(std_dh_gla.squeeze(), gsd=pleia_ddem.res[0], area=area_emp, n_patches=1000, random_state=42,
#                                                   perc_min_valid=80.)
#
#     # Second, estimate the dispersion of the means of each patch, i.e. the standard error of the mean
#     if len(df_patches)>30:
#         stderr_empirical = xdem.spatialstats.nmad(df_patches['nanmedian'].values)
#     else:
#         stderr_empirical = np.nan
#     if len(df_patches_gla) > 30:
#         stderr_empirical_gla = xdem.spatialstats.nmad(df_patches_gla['nanmedian'].values)
#     else:
#         stderr_empirical_gla = np.nan
#
#     list_stderr_empirical.append(stderr_empirical)
#     list_stderr_empirical_gla.append(stderr_empirical_gla)
#
# # Save patches results to file
# df_all_patches = pd.DataFrame()
# df_all_patches = df_all_patches.assign(area=areas_emp, stderr_emp_sta = list_stderr_empirical, stderr_emp_gla = list_stderr_empirical_gla)
# df_all_patches.to_csv('/data/icesat/travail_en_cours/romain/dem_error_study/final_run/df_patches_sta_gla.csv')

# 3/ Estimate variograms of standardized differences on stable and glacier terrain
print('FIRST VARIO')
df_vgm_sta = xdem.spatialstats.sample_empirical_variogram(pleia_ddem.data.data, pleia_ddem.res[0],
                                                          subsample=50, n_variograms=100, runs=20,
                                                          estimator='dowd', n_jobs=5, random_state=42, verbose=True)
df_vgm_sta.to_csv('/home/atom/ongoing/work_stderr_dem/case_study_montblanc/df_vgm_sta.csv')
print('SECOND VARIO')
df_vgm_gla = xdem.spatialstats.sample_empirical_variogram(pleia_ddem_gla.data.data, pleia_ddem.res[0],
                                                          subsample=50, n_variograms=100, runs=20,
                                                          estimator='dowd', n_jobs=5, random_state=42, verbose=True)
df_vgm_gla.to_csv('/home/atom/ongoing/work_stderr_dem/case_study_montblanc/df_vgm_gla.csv')
print('THIRD VARIO')
df_std_sta = xdem.spatialstats.sample_empirical_variogram(std_dh, pleia_ddem.res[0],
                                                          subsample=50, n_variograms=100, runs=20,
                                                          estimator='dowd', n_jobs=5, random_state=42, verbose=True)
df_std_sta.to_csv('/home/atom/ongoing/work_stderr_dem/case_study_montblanc/df_vgm_std_sta.csv')
print('FOURTH VARIO')
df_std_gla = xdem.spatialstats.sample_empirical_variogram(std_dh_gla, pleia_ddem.res[0],
                                                          subsample=50, n_variograms=100, runs=20,
                                                          estimator='dowd', n_jobs=5, random_state=42, verbose=True)
df_std_gla.to_csv('/home/atom/ongoing/work_stderr_dem/case_study_montblanc/df_vgm_std_gla.csv')

# Save variograms to files

# On server
# df_vgm_sta.to_csv('/data/icesat/travail_en_cours/romain/dem_error_study/final_run/df_vgm_sta.csv')
# df_vgm_gla.to_csv('/data/icesat/travail_en_cours/romain/dem_error_study/final_run/df_vgm_gla.csv')
# df_std_sta.to_csv('/data/icesat/travail_en_cours/romain/dem_error_study/final_run/df_vgm_std_sta.csv')
# df_std_gla.to_csv('/data/icesat/travail_en_cours/romain/dem_error_study/final_run/df_vgm_std_gla.csv')
