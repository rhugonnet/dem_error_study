"""Plotting of Figure S16: theoretical approximation for variogram integration"""
from typing import Callable
import matplotlib.pyplot as plt
import numpy as np
from geoutils import Raster, Vector
import xdem
from scipy.spatial.distance import pdist
import skgstat
import pandas as pd
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature

# Open data
fn_dem = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/Pleiades_Mont-Blanc_2017-10-25_DEM_5m.tif'
r = Raster(fn_dem).reproject(dst_res=200)

# Shapes with area equal to that of the Mer de Glace
fn_shp_disk = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/disk_mdg_area.shp'
fn_shp_rectangle = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/rectangle_mdg_area.shp'
fn_shp_mdg = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/mdg.shp'

fn_shp_rgi_reg11 = '/home/atom/data/inventory_products/RGI/00_rgi60_neighb_renamed/11_rgi60_CentralEurope/region_11_rgi60_CentralEurope.shp'
mdg_id = 'RGI60-11.03643'

disk = Vector(fn_shp_disk)
rectangle = Vector(fn_shp_rectangle)
rgi_reg11 = Vector(fn_shp_rgi_reg11)
mdg = Vector(rgi_reg11.ds[rgi_reg11.ds['RGIId'].isin([mdg_id])])

area_mdg = 24.179

mask_disk = disk.create_mask(r)
mask_rectangle = rectangle.create_mask(r)
mask_mdg = mdg.create_mask(r)

coords_r = np.array(r.coords())
disk_coords = coords_r[:, mask_disk.squeeze()]
rectangle_coords = coords_r[:, mask_rectangle.squeeze()]
mdg_coords = coords_r[:, mask_mdg.squeeze()]
# # radius of circle of this area
# r = np.sqrt(area/np.pi)
# # side of square of this area
# l = np.sqrt(area)
# # sides of 1x9 rectangle of this area
# short_l = l/3
# long_l = l*3

list_shape = ['Disk', 'Mer de Glace', 'Rectangle']
list_coords = [disk_coords, mdg_coords, rectangle_coords]

# Exact solution
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


# Approximate solution
def double_sum_covar_quick(coords: np.ndarray, errors: list[np.ndarray],
                     vgm_funcs: list[Callable], nb_subsample=100):
    """
    Double sum of covariance for euclidean coordinates
    :param coords: Spatial support (typically, pixel) coordinates
    :param errors: Standard errors of supports
    :param vgm_funcs: Variogram function
    :param nb_subsample: Number of points used to subset the integration
    :return:
    """

    n = len(coords)

    rand_points = np.random.choice(n, size=min(nb_subsample, n), replace=False)
    pds = pdist(coords)

    var = 0
    for ind_sub in range(nb_subsample):
        for j in range(n):

            i = rand_points[ind_sub]
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
                var += errors[k][i] * errors[k][j] * (1 - vgm_funcs[k](d))

    total_area = n * nb_subsample
    se_dsc = np.sqrt(var / total_area)

    return se_dsc

list_ranges = [400*1.2**i for i in np.arange(30)]

list_df = []
for i in range(len(list_ranges)):

    corr_range = list_ranges[i]

    print('Working on correlation range: '+str(corr_range))

    def vgm_func_short(h):
        return skgstat.models.spherical(h, corr_range, 1)

    # For MDG, rectangle and disk
    neff_rolstad = xdem.spatialstats.neff_circ(area=area_mdg * 1000000, list_vgm=[(corr_range, 'Sph', 1)])
    err_rolstad = 1/np.sqrt(neff_rolstad)

    list_err_true, list_err_approx = ([] for i in range(2))

    for c in list_coords:

        print('Working on full double sum...')
        err_true = double_sum_covar(coords=c.T, areas=np.ones(len(c.T)), errors=[np.ones(len(c.T))], vgm_funcs=[vgm_func_short])

        print('Working on approximate double sum...')
        err_approx = double_sum_covar_quick(coords=c.T, errors=[np.ones(len(c.T))], vgm_funcs=[vgm_func_short], nb_subsample=100)

        list_err_true.append(err_true)
        list_err_approx.append(err_approx)

    df_tmp = pd.DataFrame()
    df_tmp = df_tmp.assign(range=[corr_range], err_rolstad=[err_rolstad], err_true_disk=[list_err_true[0]], err_true_mdg=[list_err_true[1]],
                   err_true_rect=[list_err_true[2]], err_approx_disk=[list_err_approx[0]], err_approx_mdg=[list_err_approx[1]],
                   err_approx_rect=[list_err_approx[2]])
    list_df.append(df_tmp)

df = pd.concat(list_df)


# Initiate figure
fig = plt.figure(figsize=(10, 8))

# First, an horizontal axis on top to plot the sample histograms
ax = plt.gca()

ax.scatter(df['range'], df.err_rolstad, marker='x', color='black')
ax.scatter(df['range'], df.err_true_disk, marker='o', color='tab:orange')
ax.scatter(df['range'], df.err_true_rect, marker='o', color='tab:blue')
ax.scatter(df['range'], df.err_true_mdg, marker='o', color='tab:olive')
ax.scatter(df['range'], df.err_approx_disk, marker='<', color='tab:orange')
ax.scatter(df['range'], df.err_approx_rect, marker='<', color='tab:blue')
ax.scatter(df['range'], df.err_approx_mdg, marker='<', color='tab:olive')

ax.scatter([], [], marker='x', color='black', label='Approx. of Rolstad et al. (2009)')
ax.scatter([], [], marker='o', color='black', label='Exact integration')
ax.scatter([], [], marker='<', color='black', label='Approx. of this study')

p0 =ax.plot([], [], color='tab:orange', label='Disk shape')
p1 =ax.plot([], [], color='tab:blue', label='Rectangular shape')
p2 = ax.plot([], [], color='tab:olive', label='Mer de Glace shape')
ax.legend(loc='lower right', ncol=2)

# ax.set_yscale('logit')
ax.set_xscale('log')
ax.set_xlabel('Correlation range of variogram (m)')
ax.set_ylabel('Standardized uncertainty\nin the spatial average')

ax.text(0.5, 0.95, 'Shapes with the same area (24 km²)', transform=ax.transAxes, ha='center', va='top', fontweight='bold', fontsize=12)

ax = fig.add_axes([0, 0.4, 0.6, 0.6], projection=ccrs.UTM(32))

ax.set_extent((r.bounds.left, r.bounds.right, r.bounds.bottom, r.bounds.top), crs=ccrs.UTM(32))
ax.spines['geo'].set_visible(False)
ax.patch.set_visible(False)

shape_feature = ShapelyFeature(Reader(fn_shp_disk).geometries(), ccrs.UTM(32), edgecolor='tab:orange', alpha=1,
                               facecolor='None', linewidth=2, zorder=30)
ax.add_feature(shape_feature)


ax = fig.add_axes([0.4, 0.15, 0.6, 0.6], projection=ccrs.UTM(32))

ax.set_extent((r.bounds.left, r.bounds.right, r.bounds.bottom, r.bounds.top), crs=ccrs.UTM(32))
ax.spines['geo'].set_visible(False)
ax.patch.set_visible(False)

shape_feature = ShapelyFeature(Reader(fn_shp_rectangle).geometries(), ccrs.UTM(32), edgecolor='tab:blue', alpha=1,
                               facecolor='None', linewidth=2, zorder=30)
ax.add_feature(shape_feature)


ax = fig.add_axes([0.4, 0.3, 0.6, 0.6], projection=ccrs.UTM(32))

ax.set_extent((r.bounds.left, r.bounds.right, r.bounds.bottom, r.bounds.top), crs=ccrs.UTM(32))
ax.spines['geo'].set_visible(False)
ax.patch.set_visible(False)
shape_feature = ShapelyFeature(Reader(fn_shp_mdg).geometries(), ccrs.PlateCarree(), edgecolor='tab:olive', alpha=1,
                               facecolor='None', linewidth=2, zorder=30)
ax.add_feature(shape_feature)

# Save to file
plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_S16_final.png', dpi=400)