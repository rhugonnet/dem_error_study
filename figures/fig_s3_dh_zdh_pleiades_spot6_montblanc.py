"""Plotting of Figure S3: SPOT6-Pléiades elevation difference and standard score for the Mont-Blanc case study"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as colors
import pandas as pd
import numpy as np
import cartopy.crs as ccrs
import geoutils as gu
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
import xdem

fn_shp = '/home/atom/data/inventory_products/RGI/00_rgi60_neighb_merged/11_rgi60_CentralEurope/11_rgi60_CentralEurope.shp'
fn_df_sub = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/df_heteroscedas_slope_curv.csv'
fn_dh = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/dh_Pleiades-SPOT6_Mont-Blanc_NK_Deramp.tif'
fn_forest_shp_simplified='/home/atom/ongoing/work_stderr_dem/case_study_montblanc/outlines/forest_Mont-Blanc_ESACCI_delainey.shp'
fn_pleiades = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/Pleiades_Mont-Blanc_2017-10-25_DEM_5m.tif'

crop_ext = [339180, 5086760, 340680, 5087460]

# 1/ First, plot the elevation difference with transparent forest/glacier outlines on the full map
fig = plt.figure(figsize=(6, 5))

ax = fig.add_axes([0.025,0.375,0.45,0.6],
                      projection=ccrs.UTM(32), label='Mont-Blanc')

dh = gu.Raster(fn_dh)
plt_extent=[dh.bounds.left, dh.bounds.right, dh.bounds.bottom, dh.bounds.top]
dh_arr = gu.spatial_tools.get_array_and_mask(dh)[0]

cmap = plt.get_cmap('RdYlBu')
cmap.set_bad(color='None')

shape_feature = ShapelyFeature(Reader(fn_shp).geometries(), ccrs.PlateCarree(), edgecolor='None', alpha=0.2,
                               facecolor='tab:grey', linewidth=1, zorder=4)
ax.add_feature(shape_feature)
shape_feature = ShapelyFeature(Reader(fn_forest_shp_simplified).geometries(), ccrs.UTM(32), edgecolor='None', alpha=0.2,
                               facecolor='tab:green', linewidth=1, zorder=3)
ax.add_feature(shape_feature)

ax.imshow(dh_arr[:, :], extent=plt_extent, transform=ccrs.UTM(32), cmap=cmap,
      interpolation=None, zorder=2, vmin=-4, vmax=4)

# ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
ax.text(0.025, 0.975, 'a', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)

ax.add_patch(mpatches.Rectangle((crop_ext[0], crop_ext[1]), crop_ext[2] - crop_ext[0], crop_ext[3] - crop_ext[1],
                                edgecolor='black',facecolor='None',transform=ccrs.UTM(32),zorder=10,linewidth=1))
ax.text(crop_ext[0]+0.5*(crop_ext[2]-crop_ext[0]), crop_ext[1]-500, 'panels b & d', fontweight='bold', ha='center', va='top')


legendax = ax.inset_axes([0.9, -0.58, 0.2, 0.12], zorder=5)

legendax.set_xlim((-0.1, 0.6))
legendax.set_ylim((-0.1, 0.6))
legendax.set_xticks([])
legendax.set_yticks([])
legendax.spines['top'].set_visible(False)
legendax.spines['left'].set_visible(False)
legendax.spines['right'].set_visible(False)
legendax.spines['bottom'].set_visible(False)

legendax.add_patch(mpatches.Rectangle((0, 0), 0.2 ,0.2 , edgecolor='black',facecolor='tab:grey', alpha=0.3, zorder=10, linewidth=0.5))
legendax.text(0.25, 0.1, 'Glaciers', va='center', ha='left')
legendax.add_patch(mpatches.Rectangle((0, 0.3), 0.2 ,0.2 , edgecolor='black',facecolor='tab:green', alpha=0.3, zorder=10, linewidth=0.5))
legendax.text(0.25, 0.4, 'Forests', va='center', ha='left')

# 2/ Then, a subset zoomed on the Mer de Glace tongue
ax = fig.add_axes([0,0.12,0.5,0.225],
                      projection=ccrs.UTM(32), label='Mont-Blanc_sub')

dh.crop(crop_ext)
plt_extent=[dh.bounds.left, dh.bounds.right, dh.bounds.bottom, dh.bounds.top]
dh_arr = gu.spatial_tools.get_array_and_mask(dh)[0]
shape_feature = ShapelyFeature(Reader(fn_shp).geometries(), ccrs.PlateCarree(), edgecolor='None', alpha=0.2,
                               facecolor='tab:grey', linewidth=1, zorder=4)
ax.add_feature(shape_feature)

ax.imshow(dh_arr[:, :], extent=plt_extent, transform=ccrs.UTM(32), cmap=cmap,
      interpolation=None, zorder=2, vmin=-4, vmax=4)
ax.text(0.025, 0.975, 'b', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)


cbaxes = ax.inset_axes([0.2, -0.1, 0.6, 0.06], zorder=10)

norm = colors.Normalize(vmin=-4, vmax=4)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, cax=cbaxes, ticks=[-4, -2, 0, 2, 4], orientation='horizontal', extend='both', shrink=0.2)
cb.set_label('Elevation difference (m)')

# 3/ We standardize and replot the full map
ax = fig.add_axes([0.525,0.375,0.45,0.6],
                      projection=ccrs.UTM(32), label='Mont-Blanc')

df_sub = pd.read_csv(fn_df_sub)
fn = xdem.spatialstats.interp_nd_binning(df_sub, list_var_names=['slope_mid', 'maxc_mid'], statistic='nmad', min_count=30)

ref_dem = gu.Raster(fn_pleiades)

slope, planc, profc = xdem.terrain.get_terrain_attribute(ref_dem, attribute=['slope', 'planform_curvature',
                                                                                'profile_curvature'])
maxabsc = np.maximum(np.abs(planc), np.abs(profc))
slope_arr = gu.spatial_tools.get_array_and_mask(slope)[0]

maxabsc[maxabsc>50] = 50
dh_err = fn((slope_arr, maxabsc[0, :, :]))
dh = gu.Raster(fn_dh)
plt_extent=[dh.bounds.left, dh.bounds.right, dh.bounds.bottom, dh.bounds.top]

dh_arr = gu.spatial_tools.get_array_and_mask(dh)[0]
std_dh = dh_arr/dh_err
# std_dh[np.abs(std_dh)>7*xdem.spatialstats.nmad(std_dh)] = np.nan
std_dh /= xdem.spatialstats.nmad(std_dh)
cmap = plt.get_cmap('RdYlBu')
cmap.set_bad(color='None')

shape_feature = ShapelyFeature(Reader(fn_shp).geometries(), ccrs.PlateCarree(), edgecolor='None', alpha=0.2,
                               facecolor='tab:grey', linewidth=1, zorder=4)
ax.add_feature(shape_feature)
shape_feature = ShapelyFeature(Reader(fn_forest_shp_simplified).geometries(), ccrs.UTM(32), edgecolor='None', alpha=0.2,
                               facecolor='tab:green', linewidth=1, zorder=3)
ax.add_feature(shape_feature)

ax.imshow(std_dh[:, :], extent=plt_extent, transform=ccrs.UTM(32), cmap=cmap,
      interpolation=None, zorder=2, vmin=-3, vmax=3)

# ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
ax.text(0.025, 0.975, 'c', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)

# 4/ And plot the standardized elevation difference for the subset also
ax = fig.add_axes([0.5,0.12,0.5,0.225],
                      projection=ccrs.UTM(32), label='Mont-Blanc_sub')

std_dh_r = dh.copy(new_array=std_dh[None, :, :])
std_dh_r.crop(crop_ext)
plt_extent=[std_dh_r.bounds.left, std_dh_r.bounds.right, std_dh_r.bounds.bottom, std_dh_r.bounds.top]

std_dh_arr = gu.spatial_tools.get_array_and_mask(std_dh_r)[0]
shape_feature = ShapelyFeature(Reader(fn_shp).geometries(), ccrs.PlateCarree(), edgecolor='None', alpha=0.2,
                               facecolor='tab:grey', linewidth=1, zorder=4)
ax.add_feature(shape_feature)

ax.imshow(std_dh_arr[:, :], extent=plt_extent, transform=ccrs.UTM(32), cmap=cmap,
      interpolation=None, zorder=2, vmin=-3, vmax=3)
ax.text(0.025, 0.975, 'd', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)

cbaxes = ax.inset_axes([0.2, -0.1, 0.6, 0.06], zorder=10)

norm = colors.Normalize(vmin=-3, vmax=3)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, cax=cbaxes, ticks=[-3, -1.5, 0, 1.5, 3], orientation='horizontal', extend='both', shrink=0.2)
cb.set_label('Standardized elevation difference')

# Save to file
plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_S3_final.png', dpi=300)