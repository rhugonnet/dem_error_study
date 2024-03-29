"""Plotting of Figure 1: patterns of random and systematic errors in DEMs"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as colors
import numpy as np
import cartopy.crs as ccrs
import geoutils as gu
import xdem

# Showing dh before and after alignment with hillshade
fn_hs = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/Mont-Blanc_2017-10-25_DEM_5m_hillshade.tif'
fn_dh_nk = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/dh_shift_nk_Pleiades.tif'
fn_dh_final_dh = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/dh_Pleiades-SPOT6_Mont-Blanc_NK_Deramp.tif'

# Showing patterns of noise
fn_dh_pleiades_noise = '/home/atom/ongoing/work_stderr_dem/noise_examples/artefact_dh/dh_Peru_2017-09-01_PHR_vs_2017-08-20_PHR.tif'
fn_dh_kh9_noise = '/home/atom/ongoing/work_stderr_dem/noise_examples/artefact_dh/DZB1212-500129_003_004-DEM_coreg-diff_utm.tif'

crop_ext = [338680, 5086760, 340680, 5087460]

fig = plt.figure(figsize=(6, 5.5))

# 1/ Plot the hillshade of panel a

ax = fig.add_axes([0.25,0.79,0.5,0.2],
                      projection=ccrs.UTM(32), label='Hillshade')

hs = gu.Raster(fn_hs)
hs.crop(crop_ext)
plt_extent=[hs.bounds.left, hs.bounds.right, hs.bounds.bottom, hs.bounds.top]

color1 = colors.to_rgba('black')
color2 = colors.to_rgba('white')
cmap_ll = colors.LinearSegmentedColormap.from_list('my_cmap_hs', [color1, color2], 256)
cmap_ll.set_bad(color='None')

ax.imshow(hs.data[0, :, :], extent=plt_extent, transform=ccrs.UTM(32), cmap=cmap_ll,
      interpolation=None, zorder=2)
# ax.gridlines(draw_labels=False, dms=False, x_inline=False, y_inline=False)
ax.text(-0.1, 0.5, 'a', transform=ax.transAxes, ha='left', va='center', fontweight='bold', fontsize=14)

y_extent = hs.bounds.top - hs.bounds.bottom
x_extent = hs.bounds.right - hs.bounds.left
ax.add_patch(mpatches.Rectangle((crop_ext[2] - x_extent/20 - 400, crop_ext[1] + y_extent/5),200, 30,
                                edgecolor='black',facecolor='black',transform=ccrs.UTM(32),zorder=10,linewidth=0.5))
ax.add_patch(mpatches.Rectangle((crop_ext[2] - x_extent/20 - 200, crop_ext[1] + y_extent/5),200, 30,
                                edgecolor='black',facecolor='white',transform=ccrs.UTM(32),zorder=10,linewidth=0.5))
ax.text(crop_ext[2] - x_extent/20 - 400,  crop_ext[1] + y_extent/5 - 10,'0',ha='center',va='top',transform=ccrs.UTM(32),zorder=10)
ax.text(crop_ext[2] - x_extent/20 - 200, crop_ext[1] + y_extent/5 - 10,'200',ha='center',va='top',transform=ccrs.UTM(32),zorder=10)
ax.text(crop_ext[2] - x_extent/20 - 200, crop_ext[1] + y_extent/5 - 70,'m',ha='center',va='top',transform=ccrs.UTM(32),zorder=10)
ax.text(crop_ext[2] - x_extent/20 - 0, crop_ext[1] + y_extent/5 - 10,'400',ha='center',va='top',transform=ccrs.UTM(32),zorder=10)


# 2/ Plot the horizontal shift biases of panel a

ax = fig.add_axes([0,0.55,0.5,0.2],
                      projection=ccrs.UTM(32), label='Biases')


dh = gu.Raster(fn_dh_nk)
dh.crop(crop_ext)
hs_arr, _ = xdem.spatialstats.get_array_and_mask(dh)

cmap = plt.get_cmap('RdYlBu').copy()
cmap.set_bad(color='None')

ax.imshow(hs_arr[:, :], extent=plt_extent, transform=ccrs.UTM(32), cmap=cmap,
      interpolation=None, zorder=2, vmin=-21, vmax=-1)
ax.text(0.5, 1.02, 'Horizontal shift biases', transform=ax.transAxes, ha='center', va='bottom', fontweight='bold')

# ax.gridlines(draw_labels=False, dms=False, x_inline=False, y_inline=False)

# 3/ Plot the residuals after correcting horizontal shift

ax = fig.add_axes([0.5,0.55,0.5,0.2],
                      projection=ccrs.UTM(32), label='Heterosc')


dh_final = gu.Raster(fn_dh_final_dh)
dh_final.crop(crop_ext)
plt_extent=[dh_final.bounds.left, dh_final.bounds.right, dh_final.bounds.bottom, dh_final.bounds.top]

hs_arr, _ = xdem.spatialstats.get_array_and_mask(dh_final)
cmap = plt.get_cmap('RdYlBu').copy()
cmap.set_bad(color='None')
ax.imshow(hs_arr[:, :], extent=plt_extent, transform=ccrs.UTM(32), cmap=cmap,
      interpolation=None, zorder=2, vmin=-10, vmax=10)

# ax.gridlines(draw_labels=False, dms=False, x_inline=False, y_inline=False)
ax.text(0.5, 1.02, 'After alignment', transform=ax.transAxes, ha='center', va='bottom', fontweight='bold')


cbaxes = ax.inset_axes([-0.25, -0.125, 0.5, 0.075], zorder=10)

norm = colors.Normalize(vmin=-10, vmax=10)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, cax=cbaxes, ticks=[-10, -5, 0, 5, 10], orientation='horizontal', extend='both', shrink=0.2)
cb.set_label('Elevation difference (m)')


# 4/ Plot the along-track undulation noise from a Pléiades/Pléiades difference

ax = fig.add_axes([0.05,0.05,0.4,0.425],
                      projection=ccrs.UTM(19), label='Noise Pléiades')

dh_pleiades = gu.Raster(fn_dh_pleiades_noise)
crop_ext_pleiades = [297000, 8455000, 309000, 8463000]
dh_pleiades.crop(crop_ext_pleiades)

plt_extent=[dh_pleiades.bounds.left, dh_pleiades.bounds.right, dh_pleiades.bounds.bottom, dh_pleiades.bounds.top]

hs_arr, _ = xdem.spatialstats.get_array_and_mask(dh_pleiades)
hs_arr[np.abs(hs_arr)>5] = np.nan
cmap = plt.get_cmap('RdYlBu').copy()
cmap.set_bad(color='None')
ax.imshow(hs_arr[:, :], extent=plt_extent, transform=ccrs.UTM(19), cmap=cmap,
      interpolation=None, zorder=2, vmin=-1, vmax=1)

# ax.gridlines(draw_labels=False, dms=False, x_inline=False, y_inline=False)
ax.text(-0.1, 1.1, 'b', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)
ax.text(0.5, 1.02, 'Along-track undulations', transform=ax.transAxes, ha='center', va='bottom', fontweight='bold')

y_extent = dh_pleiades.bounds.top - dh_pleiades.bounds.bottom
x_extent = dh_pleiades.bounds.right - dh_pleiades.bounds.left
ax.add_patch(mpatches.Rectangle((crop_ext_pleiades[2] - x_extent/20 - 2000, crop_ext_pleiades[1] + y_extent/7),1000, 150,
                                edgecolor='black',facecolor='black',transform=ccrs.UTM(19),zorder=10,linewidth=0.5))
ax.add_patch(mpatches.Rectangle((crop_ext_pleiades[2] - x_extent/20 - 1000, crop_ext_pleiades[1] + y_extent/7),1000, 150,
                                edgecolor='black',facecolor='white',transform=ccrs.UTM(19),zorder=10,linewidth=0.5))
ax.text(crop_ext_pleiades[2] - x_extent/20 - 2000,  crop_ext_pleiades[1] + y_extent/7 - 100,'0',ha='center',va='top',transform=ccrs.UTM(19),zorder=10)
ax.text(crop_ext_pleiades[2] - x_extent/20 - 1000, crop_ext_pleiades[1] + y_extent/7 - 100,'1',ha='center',va='top',transform=ccrs.UTM(19),zorder=10)
ax.text(crop_ext_pleiades[2] - x_extent/20 - 1000, crop_ext_pleiades[1] + y_extent/7 - 600,'km',ha='center',va='top',transform=ccrs.UTM(19),zorder=10)
ax.text(crop_ext_pleiades[2] - x_extent/20 - 0, crop_ext_pleiades[1] + y_extent/7 - 100,'2',ha='center',va='top',transform=ccrs.UTM(19),zorder=10)


cbaxes = ax.inset_axes([0.2, -0.075, 0.6, 0.05], zorder=10)

norm = colors.Normalize(vmin=-1, vmax=1)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, cax=cbaxes, ticks=[-1, -0.5, 0, 0.5, 1], orientation='horizontal', extend='both', shrink=0.2)
cb.set_label('Elevation difference (m)')

# 5/ Plot the digitization artefacts from KH-9 imagery

ax = fig.add_axes([0.55,0.05,0.4,0.425],
                      projection=ccrs.UTM(6), label='KH9')

crop_ext_kh9 = [445000, 7634000, 505000, 7672000]


dh_kh9 = gu.Raster(fn_dh_kh9_noise)
dh_kh9.crop(crop_ext_kh9)
plt_extent=[dh_kh9.bounds.left, dh_kh9.bounds.right, dh_kh9.bounds.bottom, dh_kh9.bounds.top]
hs_arr, _ = xdem.spatialstats.get_array_and_mask(dh_kh9)

cmap = plt.get_cmap('RdYlBu')
cmap.set_bad(color='None')
ax.imshow(hs_arr[:, :], extent=plt_extent, transform=ccrs.UTM(6), cmap=cmap,
      interpolation=None, zorder=2, vmin=-10, vmax=10)

# ax.gridlines(draw_labels=False, dms=False, x_inline=False, y_inline=False)
ax.text(-0.1, 1.1, 'c', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)
ax.text(0.5, 1.02, 'Digitization artefacts', transform=ax.transAxes, ha='center', va='bottom', fontweight='bold')

y_extent = dh_kh9.bounds.top - dh_kh9.bounds.bottom
x_extent = dh_kh9.bounds.right - dh_kh9.bounds.left
ax.add_patch(mpatches.Rectangle((crop_ext_kh9[2] - x_extent/20 - 10000, crop_ext_kh9[1] + y_extent/7),5000, 750,
                                edgecolor='black',facecolor='black',transform=ccrs.UTM(6),zorder=10,linewidth=0.5))
ax.add_patch(mpatches.Rectangle((crop_ext_kh9[2] - x_extent/20 - 5000, crop_ext_kh9[1] + y_extent/7),5000, 750,
                                edgecolor='black',facecolor='white',transform=ccrs.UTM(6),zorder=10,linewidth=0.5))
ax.text(crop_ext_kh9[2] - x_extent/20 - 10000,  crop_ext_kh9[1] + y_extent/7 - 500,'0',ha='center',va='top',transform=ccrs.UTM(6),zorder=10)
ax.text(crop_ext_kh9[2] - x_extent/20 - 5000, crop_ext_kh9[1] + y_extent/7 - 500,'5',ha='center',va='top',transform=ccrs.UTM(6),zorder=10)
ax.text(crop_ext_kh9[2] - x_extent/20 - 5000, crop_ext_kh9[1] + y_extent/7 - 3000,'km',ha='center',va='top',transform=ccrs.UTM(6),zorder=10)
ax.text(crop_ext_kh9[2] - x_extent/20 - 0, crop_ext_kh9[1] + y_extent/7 - 500,'10',ha='center',va='top',transform=ccrs.UTM(6),zorder=10)


cbaxes = ax.inset_axes([0.2, -0.075, 0.6, 0.05], zorder=10)

norm = colors.Normalize(vmin=-10, vmax=10)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, cax=cbaxes, ticks=[-10, -5, 0, 5, 10], orientation='horizontal', extend='both', shrink=0.2)
cb.set_label('Elevation difference (m)')

# Save to file
plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_1_final.png', dpi=400)