"""Plotting of Figure S4: slope and curvature from the Pléiades DEM of the Mont-Blanc case study"""
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import cartopy.crs as ccrs
import geoutils as gu
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
import xdem

fn_ddem = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/dh_Pleiades-SPOT6_Mont-Blanc_NK_Deramp.tif'
fn_pleiades = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/Pleiades_Mont-Blanc_2017-10-25_DEM_5m.tif'
fn_shp = '/home/atom/data/inventory_products/RGI/00_rgi60_neighb_merged/11_rgi60_CentralEurope/11_rgi60_CentralEurope.shp'
fn_forest = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/outlines/forest_Mont-Blanc_ESACCI_delainey.shp'

pleia_ddem = gu.Raster(fn_ddem)
ref_dem = gu.Raster(fn_pleiades)
glaciers_outlines = gu.Vector(fn_shp)
forest_outlines = gu.Vector(fn_forest)
mask_glacier = glaciers_outlines.create_mask(pleia_ddem)
mask_forest = forest_outlines.create_mask(pleia_ddem)

# Remove forest, very large outliers
pleia_ddem.data[np.abs(pleia_ddem.data)>500] = np.nan

slope, planc, profc = xdem.terrain.get_terrain_attribute(ref_dem, attribute=['slope', 'planform_curvature',
                                                                                'profile_curvature'])
maxabsc = np.maximum(np.abs(planc), np.abs(profc))

# 1/ First, plot slope
fig = plt.figure(figsize=(6, 4.5))

ax = fig.add_axes([0,0.15,0.5,0.8],
                      projection=ccrs.UTM(32), label='Slope')

plt_extent=[slope.bounds.left, slope.bounds.right, slope.bounds.bottom, slope.bounds.top]
hs_arr, _ = xdem.spatialstats.get_array_and_mask(slope)

cmap = plt.get_cmap('Reds')
cmap.set_bad(color='None')

shape_feature = ShapelyFeature(Reader(fn_shp).geometries(), ccrs.PlateCarree(), edgecolor='black', alpha=0.85,
                               facecolor='None', linewidth=0.75, zorder=20)
ax.add_feature(shape_feature)
ax.imshow(hs_arr[:, :], extent=plt_extent, transform=ccrs.UTM(32), cmap=cmap,
      interpolation=None, zorder=2, vmin=0, vmax=90)

# ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False)
ax.text(0.05, 0.95, 'a', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)


cbaxes = ax.inset_axes([0.2, -0.05, 0.6, 0.025], zorder=10)

norm = colors.Normalize(vmin=0, vmax=90)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, cax=cbaxes, ticks=[0, 45, 90], orientation='horizontal', extend='both', shrink=0.2)
cb.set_label('Slope (degrees)')

# 2/ Then, plot maximum curvature
ax = fig.add_axes([0.5,0.15,0.5,0.8],
                      projection=ccrs.UTM(32), label='Curvature')

hs_arr, _ = xdem.spatialstats.get_array_and_mask(maxabsc)

col_bounds = np.array([0, 1, 5, 50])
cb = []
cb_val = np.linspace(0, 1, len(col_bounds))
cmap = plt.get_cmap('Purples')
for j in range(len(cb_val)):
    cb.append(cmap(cb_val[j]))
cmap_cus2 = colors.LinearSegmentedColormap.from_list('my_cb', list(
    zip((col_bounds - min(col_bounds)) / (max(col_bounds - min(col_bounds))), cb)), N=1000)

# cmap2._lut[0, -1] = 0.0  # We made transparent de 10 first levels of hillshade,
cmap_cus2.set_bad(color='None')

shape_feature = ShapelyFeature(Reader(fn_shp).geometries(), ccrs.PlateCarree(), edgecolor='black', alpha=0.85,
                               facecolor='None', linewidth=0.75, zorder=20)
ax.add_feature(shape_feature)
ax.imshow(hs_arr[:, :], extent=plt_extent, transform=ccrs.UTM(32), cmap=cmap_cus2,
      interpolation=None, zorder=2, vmin=0, vmax=50)

# ax.gridlines(draw_labels=False, dms=False, x_inline=False, y_inline=False)
ax.text(0.05, 0.95, 'b', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)

cbaxes = ax.inset_axes([0.2, -0.05, 0.6, 0.025], zorder=10)

norm = colors.Normalize(vmin=0, vmax=50)
sm = plt.cm.ScalarMappable(cmap=cmap_cus2, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, cax=cbaxes, ticks=[0, 5, 20, 50], orientation='horizontal', extend='both', shrink=0.2)
cb.set_label('Maximum curvature (10$^{2}$ m$^{-1}$)')

# Save to file
plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_S4_final.png', dpi=400)