"""Plotting of Figure S2: shift and tilt between Pléiades and SPOT-6 DEMs of the Mont-Blanc case study"""
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs
import geoutils as gu
import xdem

fn_shp = '/home/atom/data/inventory_products/RGI/00_rgi60_neighb_merged/11_rgi60_CentralEurope/11_rgi60_CentralEurope.shp'
fn_dh_nk = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/dh_shift_nk_Pleiades.tif'
fn_dh_deramp = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/dh_shift_deramp_Pleiades.tif'

# First, let's plot the shift
fig = plt.figure(figsize=(6, 4.25))

ax = fig.add_axes([0,0.15,0.5,0.8],
                      projection=ccrs.UTM(32), label='Mont-Blanc')

dh = gu.Raster(fn_dh_nk)
plt_extent=[dh.bounds.left, dh.bounds.right, dh.bounds.bottom, dh.bounds.top]
hs_arr, _ = xdem.spatialstats.get_array_and_mask(dh)

cmap = plt.get_cmap('RdYlBu')
cmap.set_bad(color='None')

ax.imshow(hs_arr[:, :], extent=plt_extent, transform=ccrs.UTM(32), cmap=cmap,
      interpolation=None, zorder=2, vmin=-21, vmax=-1)

# ax.gridlines(draw_labels=False, dms=False, x_inline=False, y_inline=False)
ax.text(0.05, 0.95, 'a', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)

cbaxes = ax.inset_axes([0.2, -0.05, 0.6, 0.025], zorder=10)

norm = colors.Normalize(vmin=-21, vmax=-1)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, cax=cbaxes, ticks=[-21, -16, -11, -6, -1], orientation='horizontal', extend='both', shrink=0.2)
cb.set_label('Elevation difference (m)')

# Then, plot the tilt
ax = fig.add_axes([0.5,0.15,0.5,0.8],
                      projection=ccrs.UTM(32), label='Mont-Blanc')

dh_deramp = gu.Raster(fn_dh_deramp)
plt_extent=[dh_deramp.bounds.left, dh_deramp.bounds.right, dh_deramp.bounds.bottom, dh_deramp.bounds.top]

hs_arr, _ = xdem.spatialstats.get_array_and_mask(dh_deramp)
cmap = plt.get_cmap('RdYlBu')
cmap.set_bad(color='None')
ax.imshow(hs_arr[:, :], extent=plt_extent, transform=ccrs.UTM(32), cmap=cmap,
      interpolation=None, zorder=2, vmin=-2, vmax=2)

# ax.gridlines(draw_labels=False, dms=False, x_inline=False, y_inline=False)
ax.text(0.05, 0.95, 'b', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)

cbaxes = ax.inset_axes([0.2, -0.05, 0.6, 0.025], zorder=10)

norm = colors.Normalize(vmin=-2, vmax=2)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, cax=cbaxes, ticks=[-2, -1, -0, 1, 2], orientation='horizontal', extend='both', shrink=0.2)
cb.set_label('Elevation difference (m)')

# Save to file
plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_S2_final.png', dpi=400)