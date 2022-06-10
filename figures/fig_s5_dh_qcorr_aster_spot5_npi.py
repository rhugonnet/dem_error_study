"""Plotting of Figure S5: ASTER-SPOT elevation difference and quality of stereo-correlation for the NPI case study"""
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import cartopy.crs as ccrs
import geoutils as gu
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
import xdem

fn_ddem = '/home/atom/ongoing/work_stderr_dem/case_study_npi/dh_ASTER-SPOT_NK_Deramp.tif'
fn_shp = '/home/atom/data/inventory_products/RGI/00_rgi60_neighb_merged/17_rgi60_SouthernAndes/17_rgi60_SouthernAndes.shp'
fn_corr='/home/atom/ongoing/work_stderr_dem/case_study_npi/AST_L1A_00303182012144228/AST_L1A_00303182012144228_CORR.tif'
fn_dem_aster = '/home/atom/ongoing/work_stderr_dem/case_study_npi/AST_L1A_00303182012144228/AST_L1A_00303182012144228_Z.tif'

# 1/ First, plot elevation differences
fig = plt.figure(figsize=(6, 5.5))

ax = fig.add_axes([0,0.15,0.5,0.8],
                      projection=ccrs.UTM(18, southern_hemisphere=True), label='Elevation_change')

dh = gu.Raster(fn_ddem)
dh.data[np.abs(dh.data) > 100] = np.nan
plt_extent=[dh.bounds.left, dh.bounds.right, dh.bounds.bottom, dh.bounds.top]
hs_arr, _ = xdem.spatialstats.get_array_and_mask(dh)

cmap = plt.get_cmap('RdYlBu')
cmap.set_bad(color='None')

shape_feature = ShapelyFeature(Reader(fn_shp).geometries(), ccrs.PlateCarree(), edgecolor='black', alpha=0.85,
                               facecolor='None', linewidth=0.75, zorder=20)
ax.add_feature(shape_feature)
ax.imshow(hs_arr[:, :], extent=plt_extent, transform=ccrs.UTM(18, southern_hemisphere=True), cmap=cmap,
      interpolation=None, zorder=2, vmin=-20, vmax=20)

# ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False)
ax.text(0.05, 0.95, 'a', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)


cbaxes = ax.inset_axes([0.2, -0.05, 0.6, 0.025], zorder=10)

norm = colors.Normalize(vmin=-20, vmax=20)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, cax=cbaxes, ticks=[-20, -0, 20], orientation='horizontal', extend='both', shrink=0.2)
cb.set_label('Elevation difference (m)')

# 2/ Then, plot the quality of stereo-correlation from MicMac
ax = fig.add_axes([0.5,0.15,0.5,0.8],
                      projection=ccrs.UTM(18, southern_hemisphere=True), label='Correlation')

dh_deramp = gu.Raster(fn_corr)
plt_extent=[dh_deramp.bounds.left, dh_deramp.bounds.right, dh_deramp.bounds.bottom, dh_deramp.bounds.top]

hs_arr, _ = xdem.spatialstats.get_array_and_mask(dh_deramp)
cmap = plt.get_cmap('Greens')
cmap.set_bad(color='None')

shape_feature = ShapelyFeature(Reader(fn_shp).geometries(), ccrs.PlateCarree(), edgecolor='black', alpha=0.85,
                               facecolor='None', linewidth=0.75, zorder=20)
ax.add_feature(shape_feature)
ax.imshow(hs_arr[:, :], extent=plt_extent, transform=ccrs.UTM(18, southern_hemisphere=True), cmap=cmap,
      interpolation=None, zorder=2, vmin=0, vmax=100)

# ax.gridlines(draw_labels=False, dms=False, x_inline=False, y_inline=False)
ax.text(0.05, 0.95, 'b', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)

cbaxes = ax.inset_axes([0.2, -0.05, 0.6, 0.025], zorder=10)

norm = colors.Normalize(vmin=0, vmax=100)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, cax=cbaxes, ticks=[0, 50, 100], orientation='horizontal', extend='both', shrink=0.2)
cb.set_label('Quality of stereo-correlation (%)')

# Save to file
plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_S5_final.png', dpi=300)