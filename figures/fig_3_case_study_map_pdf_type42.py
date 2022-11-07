"""Plotting of Figure 3: map of Mont-Blanc case study"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as colors
import os
import cartopy.crs as ccrs
import geoutils as gu
import pyproj
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
import cartopy.feature as cfeature

plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'pdf.fonttype':42})
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["text.usetex"] = True

fn_lc = '/home/atom/data/inventory_products/Land_cover/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7_robinson.tif'
fn_hs = '/home/atom/documents/paper/Hugonnet_2020/figures/world_robin_rs.tif'
fn_land = '/home/atom/data/inventory_products/NaturalEarth/ne_50m_land/ne_50m_land.shp'
fn_shp = '/home/atom/data/inventory_products/RGI/00_rgi60_neighb_merged/11_rgi60_CentralEurope/11_rgi60_CentralEurope.shp'
fn_shp_buff = '/home/atom/data/inventory_products/RGI/buffered/rgi60_buff_diss.shp'
fn_ais = '/home/atom/data/inventory_products/RGI/AIS_mask/ais_glacier_ice_mask_wgs84.shp'
fn_gis = '/home/atom/data/inventory_products/RGI/GIS_mask/GreenlandMasks/Greenland_IceMask_wgs84.shp'
fn_forest_shp_simplified='/home/atom/ongoing/work_stderr_dem/case_study_montblanc/outlines/forest_Mont-Blanc_ESACCI_delainey.shp'

fn_hs_montblanc = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/Mont-Blanc_2017-10-25_DEM_5m_hillshade.tif'

hs_r = gu.Raster(fn_hs)

fig = plt.figure(figsize=(6, 7))

grid = plt.GridSpec(20, 20, wspace=0.1, hspace=0.1)

# First panel: world map with coverage + inset for Mont-Blanc case study
# sub_ax = fig.add_axes([0,0.1,0.6,0.9],
#                       projection=ccrs.Robinson(), label='world')
#
# bounds = [-179.99,179.99,-89.99,89.99]
#
# def poly_from_extent(ext):
#
#     poly = np.array([(ext[0],ext[2]),(ext[1],ext[2]),(ext[1],ext[3]),(ext[0],ext[3]),(ext[0],ext[2])])
#
#     return poly
# polygon = poly_from_extent(bounds)
#
# # Add hillshade
# img = GeoImg(fn_hs)
# land_mask = create_mask_from_shapefile(img,fn_land)
# ds = gdal.Open(fn_hs)
# gt = ds.GetGeoTransform()  # Defining bounds
# ext = (gt[0], gt[0] + ds.RasterXSize * gt[1],
#            gt[3] + ds.RasterYSize * gt[5], gt[3])
# hs = ds.ReadAsArray()
# hs = hs.astype(float)
# ds = None
#
# def stretch_hs(hs,stretch_factor=1.):
#
#     max_hs = 255
#     min_hs = 0
#
#     hs_s = (hs - (max_hs-min_hs)/2)*stretch_factor + (max_hs-min_hs)/2
#
#     return hs_s
#
# hs = stretch_hs(hs,stretch_factor=0.9)
#
# hs_land = hs.copy()
# hs_land[~land_mask]=0
# hs_notland = hs.copy()
# hs_notland[land_mask]=0
#
# hs_tmp = hs_land.copy()
# hs_tmp_nl = hs_notland.copy()
#
# def inter_poly_coords(polygon_coords):
#     list_lat_interp = []
#     list_lon_interp = []
#     for i in range(len(polygon_coords) - 1):
#         lon_interp = np.linspace(polygon_coords[i][0], polygon_coords[i + 1][0], 50)
#         lat_interp = np.linspace(polygon_coords[i][1], polygon_coords[i + 1][1], 50)
#
#         list_lon_interp.append(lon_interp)
#         list_lat_interp.append(lat_interp)
#
#     all_lon_interp = np.concatenate(list_lon_interp)
#     all_lat_interp = np.concatenate(list_lat_interp)
#
#     return np.array(list(zip(all_lon_interp, all_lat_interp)))
#
# def out_of_poly_mask(geoimg, poly_coords):
#
#     poly = ot.poly_from_coords(inter_poly_coords(poly_coords))
#     srs = osr.SpatialReference()
#     srs.ImportFromEPSG(4326)
#
#     # put in a memory vector
#     ds_shp = ot.create_mem_shp(poly, srs)
#
#     return ot.geoimg_mask_on_feat_shp_ds(ds_shp, geoimg)
#
# mask = out_of_poly_mask(img, polygon)
#
# hs_tmp[~mask] = 0
# hs_tmp_nl[~mask] = 0
#
# color1 = mpl.colors.to_rgba('black')
# color2 = mpl.colors.to_rgba('white')
# cmap2 = mpl.colors.LinearSegmentedColormap.from_list('my_cmap2', [color1, color2], 256)
# cmap2._init()
# cmap2._lut[0:1, -1] = 0.0  # We made transparent de 10 first levels of hillshade,
# cmap2._lut[1:, -1] = 0.60
#
# cmap22 = mpl.colors.LinearSegmentedColormap.from_list('my_cmap22', [color1, color2], 256)
# cmap22._init()
# cmap22._lut[0:1, -1] = 0.0  # We made transparent de 10 first levels of hillshade,
# cmap22._lut[1:, -1] = 0.3
#
# sc_img = GeoImg(fn_out)
# cmap_sc = mpl.colors.LinearSegmentedColormap.from_list('my_cmap_sc', [color1, color2], 2)
# cmap_sc._init()
# cmap_sc._lut[0:1, -1] = 0
# cmap_sc._lut[1, -1] = 0.9
#
# lc_r = gu.Raster(fn_lc)
# lc_arr, _ = gu.spatial_tools.get_array_and_mask(lc_r)
# forest_mask = np.logical_or(np.logical_and(lc_arr>=50, lc_arr<=90), lc_arr==160, lc_arr==170)
# water_bodies = lc_arr==210
# water_bodies[~land_mask] = False
#
# cmap_forest = mpl.colors.LinearSegmentedColormap.from_list('my_cmap_forest', [ color1, mpl.colors.to_rgba('tab:green')], 2)
# cmap_forest._init()
# cmap_forest._lut[0:1, -1] = 0
# cmap_forest._lut[1, -1] = 0.8
#
# cmap_wb = mpl.colors.LinearSegmentedColormap.from_list('my_cmap_waterbodies', [ color1, mpl.colors.to_rgba('tab:blue')], 2)
# cmap_wb._init()
# cmap_wb._lut[0:1, -1] = 0
# cmap_wb._lut[1, -1] = 0.8
#
#
# shape_feature = ShapelyFeature(Reader(fn_shp_buff).geometries(), ccrs.PlateCarree(), edgecolor='None', alpha=0.5,
#                                facecolor='tab:cyan', linewidth=0, zorder=6)
# sub_ax.add_feature(shape_feature)
# shape_feature = ShapelyFeature(Reader(fn_ais).geometries(), ccrs.PlateCarree(), edgecolor='None', alpha=0.5,
#                                facecolor='tab:cyan', linewidth=0, zorder=6)
# sub_ax.add_feature(shape_feature)
# shape_feature = ShapelyFeature(Reader(fn_gis).geometries(), ccrs.PlateCarree(), edgecolor='None', alpha=0.5,
#                                facecolor='tab:cyan', linewidth=0, zorder=6)
# sub_ax.add_feature(shape_feature)
#
# shape_feature = ShapelyFeature(Reader(fn_studysites).geometries(), ccrs.PlateCarree(), edgecolor='none', alpha=0.45,
#                                facecolor='tab:orange', linewidth=0, zorder=7)
# sub_ax.add_feature(shape_feature)
# shape_feature = ShapelyFeature(Reader(fn_studysites).geometries(), ccrs.PlateCarree(), edgecolor='tab:orange', alpha=1,
#                                facecolor='none', linewidth=0.5, zorder=7)
# sub_ax.add_feature(shape_feature)
#
# sub_ax.imshow(hs_tmp[:, :], extent=ext, transform=ccrs.Robinson(), cmap=cmap2, zorder=2, interpolation='nearest',rasterized=True)
# sub_ax.imshow(hs_tmp_nl[:, :], extent=ext, transform=ccrs.Robinson(), cmap=cmap22, zorder=2,interpolation='nearest',rasterized=True)
#
# sub_ax.imshow(sc_img.img[:, :], extent=ext, transform=ccrs.Robinson(), cmap=cmap_sc, zorder=3, interpolation='nearest', rasterized=True)
#
# sub_ax.imshow(forest_mask.astype(np.float32)[:, :], extent=ext, transform=ccrs.Robinson(), cmap=cmap_forest, zorder=4, interpolation='nearest', rasterized=True)
# sub_ax.imshow(water_bodies.astype(np.float32)[:, :], extent=ext, transform=ccrs.Robinson(), cmap=cmap_wb, zorder=5, interpolation='nearest', rasterized=True)
#
#
# sub_ax.set_extent([-179.99,179.99,-89.99,89.99], ccrs.Geodetic())
#
# sub_ax.add_feature(cfeature.NaturalEarthFeature('physical', 'ocean', '50m', facecolor='lightgrey'), alpha=0.5)
# sub_ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', facecolor=plt.cm.Greys(0.9)), alpha=0.5)
#
# sub_ax.outline_patch.set_edgecolor('lightgrey')
# sub_ax.text(0.05, 0.95, 'a', transform=sub_ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)


# Add the coverage of various data
# shape_feature = ShapelyFeature(Reader(fn_shp).geometries(), ccrs.PlateCarree(), edgecolor='None', alpha=0.5,
#                                facecolor='tab:cyan', linewidth=1)
# sub_ax.add_feature(shape_feature)


# 2/ Legend (plot in advance to be behind)

legendax = fig.add_axes([0, 0.035, 0.2, 0.2], label='legends')

legendax.set_xlim((-0.5, 0.5))
legendax.set_ylim((-0.15, 1.15))
legendax.set_xticks([])
legendax.set_yticks([])
legendax.spines['top'].set_visible(False)
legendax.spines['left'].set_visible(False)
legendax.spines['right'].set_visible(False)
legendax.spines['bottom'].set_visible(False)

legendax.add_patch(mpatches.Rectangle((0, 0.25), 0.2 ,0.2 , edgecolor='black',facecolor='tab:cyan', zorder=10, linewidth=0.5))
legendax.text(0.3, 0.35, 'Glacierized', va='center', ha='left')

legendax = fig.add_axes([0.575, 0.035, 0.2, 0.2], label='legends3')

legendax.set_xlim((-0.5, 0.5))
legendax.set_ylim((-0.15, 1.15))
legendax.set_xticks([])
legendax.set_yticks([])
legendax.spines['top'].set_visible(False)
legendax.spines['left'].set_visible(False)
legendax.spines['right'].set_visible(False)
legendax.spines['bottom'].set_visible(False)

legendax.add_patch(mpatches.Rectangle((0, 0.25), 0.2 ,0.2 , edgecolor='darkblue',facecolor='white', zorder=10, linewidth=2))
legendax.text(0.3, 0.35, 'Example glaciers', va='center', ha='left')

legendax = fig.add_axes([0.3, 0.035, 0.2, 0.2], label='legends2')

legendax.set_xlim((-0.5, 0.5))
legendax.set_ylim((-0.15, 1.15))
legendax.set_xticks([])
legendax.set_yticks([])
legendax.spines['top'].set_visible(False)
legendax.spines['left'].set_visible(False)
legendax.spines['right'].set_visible(False)
legendax.spines['bottom'].set_visible(False)

legendax.add_patch(mpatches.Rectangle((0, 0.25), 0.2 ,0.2 , edgecolor='black',facecolor='tab:green', zorder=10, linewidth=0.5))
legendax.text(0.3, 0.35, 'Forested', va='center', ha='left')

# 3/ Plot Mont-Blanc case study map

ax = fig.add_axes([0.075,0.15,0.9,0.8],
                      projection=ccrs.UTM(32), label='Mont-Blanc')

hs = gu.Raster(fn_hs_montblanc)
plt_extent=[hs.bounds.left, hs.bounds.right, hs.bounds.bottom, hs.bounds.top]
crop_extent = [plt_extent[0], plt_extent[2], plt_extent[1], plt_extent[3]]

hs_arr = hs.data

ax.set_extent([hs.bounds.left -1500, hs.bounds.right+1500, hs.bounds.bottom, hs.bounds.top], ccrs.UTM(32))

color1 = colors.to_rgba('black')
color2 = colors.to_rgba('white')
cmap_ll = colors.LinearSegmentedColormap.from_list('my_cmap_hs', [color1, color2], 256)
# cmap_ll._init()
# cmap_ll._lut[1:, -1] = 1
# cmap_ll._lut[0:1, -1] = 0.0  # We make transparent the lowest hillshade value
cmap_ll.set_bad(color='None')

shape_feature = ShapelyFeature(Reader(fn_shp).geometries(), ccrs.PlateCarree(), edgecolor='None', alpha=0.65,
                               facecolor='tab:cyan', linewidth=1, zorder=4)
ax.add_feature(shape_feature)

ax.imshow(hs_arr[0, :, :], extent=plt_extent, transform=ccrs.UTM(32), cmap=cmap_ll,
      interpolation=None, zorder=2)

y_extent = hs.bounds.top - hs.bounds.bottom
x_extent = hs.bounds.right - hs.bounds.left
ax.add_patch(mpatches.Rectangle((crop_extent[2] - x_extent/20 - 4000, crop_extent[1] + y_extent/15),2000, 300,
                                edgecolor='black',facecolor='black',transform=ccrs.UTM(32),zorder=10,linewidth=0.5))
ax.add_patch(mpatches.Rectangle((crop_extent[2] - x_extent/20 - 2000, crop_extent[1] + y_extent/15),2000, 300,
                                edgecolor='black',facecolor='white',transform=ccrs.UTM(32),zorder=10,linewidth=0.5))
ax.text(crop_extent[2] - x_extent/20 - 4000,  crop_extent[1] + y_extent/15 - 100,'0',ha='center',va='top',transform=ccrs.UTM(32),zorder=10)
ax.text(crop_extent[2] - x_extent/20 - 2000, crop_extent[1] + y_extent/15 - 100,'2',ha='center',va='top',transform=ccrs.UTM(32),zorder=10)
ax.text(crop_extent[2] - x_extent/20 - 2000, crop_extent[1] + y_extent/15 - 800,'km',ha='center',va='top',transform=ccrs.UTM(32),zorder=10)
ax.text(crop_extent[2] - x_extent/20 - 0, crop_extent[1] + y_extent/15 - 100,'4',ha='center',va='top',transform=ccrs.UTM(32),zorder=10)

shape_feature = ShapelyFeature(Reader(fn_forest_shp_simplified).geometries(), ccrs.UTM(32), edgecolor='None', alpha=0.65,
                               facecolor='tab:green', linewidth=1, zorder=3)
ax.add_feature(shape_feature)

ax.gridlines(draw_labels={'top':'x', 'left':'y'}, dms=True, x_inline=False, y_inline=False)

ll = (hs.bounds.left, hs.bounds.bottom)
lr = (hs.bounds.right, hs.bounds.bottom)
ul = (hs.bounds.left, hs.bounds.top)
ur = (hs.bounds.right, hs.bounds.top)

# Extent of Figure on heteroscedasticty
y_extent = hs.bounds.top - hs.bounds.bottom
x_extent = hs.bounds.right - hs.bounds.left
plt_extent0 = [
    hs.bounds.left + 3/10*x_extent,
    hs.bounds.right - 2/10*x_extent,
    hs.bounds.bottom + 1/10*y_extent,
    hs.bounds.top - 4.3/10*y_extent,
]

ax.add_patch(mpatches.Rectangle((plt_extent0[0], plt_extent0[2]), plt_extent0[1] - plt_extent0[0], plt_extent0[3] - plt_extent0[2],
                                edgecolor='black',facecolor='none',transform=ccrs.UTM(32),zorder=10,linewidth=2, linestyle='dashed'))

ax.text(plt_extent0[0]+500, plt_extent0[3]-500,
        'Fig. 4(d)', ha='left', va='top', bbox=dict(edgecolor='black', facecolor='white', alpha=1), zorder=5)

# Extent of Figure on spatial correlations
plt_extent = [
    hs.bounds.left+1/40*x_extent,
    hs.bounds.right-1/5*x_extent,
    hs.bounds.bottom + 6/10*y_extent,
    hs.bounds.top-1/20*y_extent,
]
ax.add_patch(mpatches.Rectangle((plt_extent[0], plt_extent[2]), plt_extent[1] - plt_extent[0], plt_extent[3] - plt_extent[2],
                                edgecolor='black',facecolor='none',transform=ccrs.UTM(32),zorder=10,linewidth=2, linestyle='dashed'))
ax.text(plt_extent[1]-500, plt_extent[3]-500, 'Fig. 5(b)', ha='right', va='top', bbox=dict(edgecolor='black', facecolor='white', alpha=1), zorder=5)

# Extent of Figure on slope simulation
crop_ext = [333500, 5076000, 335500, 5078000]

ax.add_patch(mpatches.Rectangle((crop_ext[0], crop_ext[1]), crop_ext[2] - crop_ext[0], crop_ext[3] - crop_ext[1],
                                edgecolor='black',facecolor='none',transform=ccrs.UTM(32),zorder=10,linewidth=2, linestyle='dashed'))
ax.text(crop_ext[2]+500, crop_ext[3]-900, 'Fig. 6(a)', ha='left', va='top', bbox=dict(edgecolor='black', facecolor='white', alpha=1), zorder=5)

# Display glaciers used in Table 2
names_gla = ['Bossons', 'Taconnaz', 'Griaz', 'Bourgeat']
rgiid_gla = ['RGI60-11.03646', 'RGI60-11.03647', 'RGI60-11.03280', 'RGI60-11.03290']
glacier_inventory = gu.Vector(fn_shp)
shift_yy = [3100, 2300, -100, 900]
shift_xx = [-1500, -900, -1700, -800]
for rgiid in rgiid_gla:
    gla_ds = glacier_inventory.ds[glacier_inventory.ds['RGIId'].values==rgiid]
    shape_feature = ShapelyFeature(gla_ds.geometry, ccrs.PlateCarree(), edgecolor='darkblue',
                                   alpha=1,
                                   facecolor='None', linewidth=1, zorder=4)
    ax.add_feature(shape_feature)
    lat, lon = gla_ds['CenLat'].values[0], gla_ds['CenLon'].values[0]
    xx, yy = gu.projtools.reproject_from_latlon((lat, lon), out_crs=pyproj.CRS(32632))

    ax.text(xx + shift_xx[rgiid_gla.index(rgiid)], yy + shift_yy[rgiid_gla.index(rgiid)], '$$\\textbf{'+names_gla[rgiid_gla.index(rgiid)]+'}$$', va='center', ha='center', color='darkblue', zorder=30, fontweight='bold')


# Mont-Blanc location
x_mb = 334360
y_mb = 5077742
ax.scatter(x_mb, y_mb, s=100, marker='x', color='white', linewidths=3, zorder=30)
ax.text(x_mb - 800 , y_mb, '$$\\textbf{Mont}$$'+'\n'+'$$\\textbf{Blanc}$$', color='white', va='center', ha='right', zorder=30)

# ax.text(0.025, 0.975, 'b', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14, zorder=30)

# 4/ Add inset

sub_ax = fig.add_axes([0.05,0.775,0.25,0.25],
                      projection=ccrs.Robinson(), label='world')

sub_ax.set_facecolor("white")

sub_ax.add_feature(cfeature.NaturalEarthFeature('physical', 'ocean', '50m', facecolor='white'), alpha=1)
sub_ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', facecolor='white'), alpha=1)

sub_ax.add_feature(cfeature.NaturalEarthFeature('physical', 'ocean', '50m', facecolor='lightgrey'), alpha=0.3)
sub_ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', facecolor=plt.cm.Greys(0.95)), alpha=0.8)
plt_extent=[hs_r.bounds.left, hs_r.bounds.right, hs_r.bounds.bottom, hs_r.bounds.top]

hs_arr = gu.spatial_tools.get_array_and_mask(hs_r)[0]
color1 = colors.to_rgba('black')
color2 = colors.to_rgba('white')
cmap2 = colors.LinearSegmentedColormap.from_list('my_cmap2', [color1, color2], 256)
cmap2._init()
cmap2._lut[0:1, -1] = 0.0  # We made transparent de 10 first levels of hillshade,
cmap2._lut[1:, -1] = 0.60
sub_ax.imshow(hs_arr[:, :], transform=ccrs.Robinson(), extent=plt_extent, cmap=cmap2, zorder=2, interpolation='nearest',rasterized=True)

bounds = [-10, 20, 35, 55]
sub_ax.set_extent(bounds, ccrs.Geodetic())

sub_ax.plot(6.86, 45.83, marker='s', color='red', transform=ccrs.Geodetic())
# sub_ax.text(10, 50, 'Europe', color=plt.cm.Greys(0.8), ha='center', va='center', transform=ccrs.Geodetic(), fontweight='bold', rotation=15)

# Save to file
plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_3_final.pdf', dpi=400, transparent=True)