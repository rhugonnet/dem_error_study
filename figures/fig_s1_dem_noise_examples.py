"""Plotting of Figure S1: example of correlated noises in DEMs"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as colors
import cartopy.crs as ccrs
import geoutils as gu
import xdem

fn_shp = '/home/atom/data/inventory_products/RGI/00_rgi60_neighb_merged/11_rgi60_CentralEurope/11_rgi60_CentralEurope.shp'

# Showing patterns of noise
fn_dh_spot6_glo_noise = '/home/atom/ongoing/work_stderr_dem/noise_examples/artefact_dh/dhdt_Iceland_2020-08-08_SPOT6_vs_GLO30.tif'
fn_dh_srtm_x = '/home/atom/ongoing/work_stderr_dem/noise_examples/artefact_dh/dh-HGTS-XSAR/example_noise.vrt'
fn_dh_wv_tdx = '/home/atom/ongoing/work_stderr_dem/noise_examples/dh_TDX_WV.tif'
fn_dh_ast_tdx = '/home/atom/ongoing/work_stderr_dem/noise_examples/AST_L1A_00311202000201156/dh_AST_TDX_Nahanni.tif'

# COMMENTED: how the difference were processed from original segments

# fn_wv = '/home/atom/ongoing/work_stderr_dem/final/noise_examples/SETSM_WV02_20121120_103001001DAEB200_103001001CC8DE00_seg12_2m_v3.0/SETSM_WV02_20121120_103001001DAEB200_103001001CC8DE00_seg12_2m_v3.0_dem.tif'
# fn_tdx = '/home/atom/ongoing/work_stderr_dem/final/noise_examples/TDX_90m_05hem_N61W128.tif'
# wv = xdem.DEM(fn_wv)
# tdx = xdem.DEM(fn_tdx)
# wv = wv.reproject(tdx, resampling=rio.enums.Resampling.bilinear)
# dh_wv_tdx = wv - tdx
# dh_wv_tdx.save('/home/atom/ongoing/work_stderr_dem/final/noise_examples/dh_TDX_WV.tif')
# dh_wv_tdx.show(vmin=-10, vmax=10, cmap='Spectral')

# fn_aster = '/home/atom/ongoing/work_stderr_dem/final/noise_examples/AST_L1A_00311202000201156/AST_L1A_00311202000201156_Z.tif'
# fn_tdx = '/home/atom/ongoing/work_stderr_dem/final/noise_examples/TDX_for_ASTER_example/tdx.vrt'
#
# aster = xdem.DEM(fn_aster)
# tdx = xdem.DEM(fn_tdx)
# tdx = tdx.reproject(aster, resampling=rio.enums.Resampling.bilinear)
#
# coreg = xdem.coreg.NuthKaab()
# coreg.fit(reference_dem=tdx, dem_to_be_aligned=aster, verbose=True)
# aster_aligned = coreg.apply(aster)
# dh = aster_aligned - tdx
# corr = gu.Raster('/home/atom/ongoing/work_stderr_dem/final/noise_examples/AST_L1A_00311202000201156/AST_L1A_00311202000201156_CORR.tif')
# dh.data[corr.data<60.]=np.nan
#
# dh.save('/home/atom/ongoing/work_stderr_dem/final/noise_examples/AST_L1A_00311202000201156/dh_AST_TDX_Nahanni.tif')


# Function for distance scale bar on maps
def scale_bar(r, distance, ax, distance_labels):

    utm = str(dh.crs.to_epsg())[-2:]
    s_n = str(dh.crs.to_epsg())[-3]
    crs = ccrs.UTM(utm, southern_hemisphere=s_n == '7')

    ax.add_patch(mpatches.Rectangle((r.bounds.right - distance*1.5, r.bounds.bottom + distance*0.8), distance, 0.15*distance,
                                    edgecolor='black',facecolor='black',transform=crs,zorder=10,linewidth=0.5))
    ax.add_patch(mpatches.Rectangle((r.bounds.right - distance*1.5 - distance, r.bounds.bottom + distance*0.8), distance, 0.15*distance,
                                    edgecolor='black',facecolor='white',transform=crs ,zorder=10,linewidth=0.5))
    ax.text(r.bounds.right - distance*1.5 - distance,  r.bounds.bottom + distance*0.7,distance_labels[0],ha='center',va='top',transform=crs,zorder=10)
    ax.text(r.bounds.right - distance*1.5, r.bounds.bottom + distance*0.7, distance_labels[1],ha='center',va='top',transform=crs,zorder=10)
    ax.text(r.bounds.right - distance*1.5, r.bounds.bottom + distance*0.4,distance_labels[3],ha='center',va='top',transform=crs,zorder=10)
    ax.text(r.bounds.right - distance*1.5 + distance, r.bounds.bottom + distance*0.7,distance_labels[2],ha='center',va='top',transform=crs,zorder=10)

# 1/ First, ASTER undulations
fig = plt.figure(figsize=(12, 12))

ax = fig.add_axes([0.025,0.6,0.45,0.35],
                      projection=ccrs.UTM(9), label='ASTER')


dh = gu.Raster(fn_dh_ast_tdx)
crop_ext = [525000, 6653000, 606000, 6713000]
dh.crop(crop_ext)
hs_arr, _ = xdem.spatialstats.get_array_and_mask(dh)

plt_extent = [crop_ext[0], crop_ext[2], crop_ext[1], crop_ext[3]]

cmap = plt.get_cmap('RdYlBu').copy()
cmap.set_bad(color='None')

ax.imshow(hs_arr[:, :], extent=plt_extent, transform=ccrs.UTM(9), cmap=cmap,
      interpolation=None, zorder=2, vmin=-15, vmax=15)
ax.text(0.5, 1.02, 'ASTER: cross-track bias and along-track undulations', transform=ax.transAxes, ha='center', va='bottom', fontweight='bold')
ax.text(-0.05, 1.05, 'a', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)

ax.gridlines(draw_labels=False, dms=False, x_inline=False, y_inline=False)

scale_bar(dh, 10000, ax=ax, distance_labels=['0', '10', '20', 'km'])

cbaxes = ax.inset_axes([0.2, -0.075, 0.6, 0.035], zorder=10)

norm = colors.Normalize(vmin=-15, vmax=15)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, cax=cbaxes, ticks=[-15, -0, 15], orientation='horizontal', extend='both', shrink=0.2)
cb.set_label('Elevation difference (m)')


# 2/ Second, SRTM-X undulations
ax = fig.add_axes([0.525,0.6,0.45,0.35],
                      projection=ccrs.UTM(45), label='SRTMX')


dh = gu.Raster(fn_dh_srtm_x)
crop_ext = [680000, 3765000, 1138000, 4105000]
dh.crop(crop_ext)
hs_arr, _ = xdem.spatialstats.get_array_and_mask(dh)

plt_extent = [crop_ext[0], crop_ext[2], crop_ext[1], crop_ext[3]]

cmap = plt.get_cmap('RdYlBu').copy()
cmap.set_bad(color='None')

ax.imshow(hs_arr[:, :], extent=plt_extent, transform=ccrs.UTM(45), cmap=cmap,
      interpolation=None, zorder=2, vmin=-5, vmax=5)
ax.text(0.5, 1.02, 'SRTM-X: along-track undulations', transform=ax.transAxes, ha='center', va='bottom', fontweight='bold')
ax.text(-0.05, 1.05, 'b', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)

ax.gridlines(draw_labels=False, dms=False, x_inline=False, y_inline=False)

scale_bar(dh, 50000, ax=ax, distance_labels=['0', '50', '100', 'km'])

cbaxes = ax.inset_axes([0.2, -0.075, 0.6, 0.035], zorder=10)

norm = colors.Normalize(vmin=-5, vmax=5)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, cax=cbaxes, ticks=[-5, -0, 5], orientation='horizontal', extend='both', shrink=0.2)
cb.set_label('Elevation difference (m)')

# 3/ Third, Worldview artefacts
ax = fig.add_axes([0.025,0.1,0.45,0.35],
                      projection=ccrs.UTM(9), label='WV')


dh = gu.Raster(fn_dh_wv_tdx)
crop_ext = [580000, 6861500, 583000, 6864000]
dh.crop(crop_ext)
hs_arr, _ = xdem.spatialstats.get_array_and_mask(dh)

plt_extent = [crop_ext[0], crop_ext[2], crop_ext[1], crop_ext[3]]

cmap = plt.get_cmap('RdYlBu').copy()
cmap.set_bad(color='None')

ax.imshow(hs_arr[:, :], extent=plt_extent, transform=ccrs.UTM(9), cmap=cmap,
      interpolation=None, zorder=2, vmin=-50, vmax=50)

ax.gridlines(draw_labels=False, dms=False, x_inline=False, y_inline=False)
ax.text(-0.05, 1.05, 'c', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)
ax.text(0.5, 1.02, 'ArcticDEM: processing artefacts', transform=ax.transAxes, ha='center', va='bottom', fontweight='bold')

scale_bar(dh, 250, ax=ax, distance_labels=['0', '250', '500', 'm'])

cbaxes = ax.inset_axes([0.2, -0.075, 0.6, 0.035], zorder=10)

norm = colors.Normalize(vmin=-50, vmax=50)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, cax=cbaxes, ticks=[-50, 0, 50], orientation='horizontal', extend='both', shrink=0.2)
cb.set_label('Elevation difference (m)')

# 4/ Finally, SPOT6 artefacts
ax = fig.add_axes([0.525,0.1,0.45,0.35],
                      projection=ccrs.UTM(28), label='SPOT6')

dh = gu.Raster(fn_dh_spot6_glo_noise)
dh.set_ndv(-9999)
crop_ext = [449000, 7115000, 459000, 7123000]
dh.crop(crop_ext)
hs_arr, _ = xdem.spatialstats.get_array_and_mask(dh)

plt_extent = [crop_ext[0], crop_ext[2], crop_ext[1], crop_ext[3]]

cmap = plt.get_cmap('RdYlBu').copy()
cmap.set_bad(color='None')

ax.imshow(hs_arr[:, :], extent=plt_extent, transform=ccrs.UTM(28), cmap=cmap,
      interpolation=None, zorder=2, vmin=-10, vmax=10)
ax.text(0.5, 1.02, 'SPOT-6: processing artefacts', transform=ax.transAxes, ha='center', va='bottom', fontweight='bold')
ax.text(-0.05, 1.05, 'd', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)

ax.gridlines(draw_labels=False, dms=False, x_inline=False, y_inline=False)
scale_bar(dh, 1000, ax=ax, distance_labels=['0', '1', '2', 'km'])

cbaxes = ax.inset_axes([0.2, -0.075, 0.6, 0.035], zorder=10)

norm = colors.Normalize(vmin=-10, vmax=10)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, cax=cbaxes, ticks=[-10, -0, 10], orientation='horizontal', extend='both', shrink=0.2)
cb.set_label('Elevation difference (m)')

# Save to file
plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_S1_final.png', dpi=400)