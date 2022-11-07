"""Pre-process Pléiades and SPOT-6 DEMs into elevation changes for the Mont-Blanc case study"""
import os
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    import xdem
import geoutils as gu
import numpy as np

# Open DEMs and outlines
fn_glacier_outlines = '/home/atom/data/inventory_products/RGI/00_rgi60_neighb_renamed/11_rgi60_CentralEurope/region_11_rgi60_CentralEurope.shp'
fn_dem_spot = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/SPOT6_Mont-Blanc_2017-10-24_DEM_5m.tif'
fn_dem_pleiades = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/Pleiades_Mont-Blanc_2017-10-25_DEM_5m.tif'
fn_forest_shp_simplified='/home/atom/ongoing/work_stderr_dem/case_study_montblanc/outlines/forest_Mont-Blanc_ESACCI_delainey.shp'

dem_spot = xdem.DEM(fn_dem_spot)
dem_pleiades = xdem.DEM(fn_dem_pleiades)
glacier_outlines = gu.Vector(fn_glacier_outlines)
forest_outlines = gu.Vector(fn_forest_shp_simplified)

# Reproject on the Pléiades DEM
reproj_dem_spot = dem_spot.reproject(dem_pleiades)
init_dh = reproj_dem_spot.data - dem_pleiades.data

# Create mask of inlier data
inlier_mask = np.logical_and.reduce((~glacier_outlines.create_mask(dem_pleiades),
                                     ~forest_outlines.create_mask(dem_pleiades),
                                     np.abs(init_dh - np.nanmedian(init_dh))<3*xdem.spatialstats.nmad(init_dh)))

# Coregistration pipeline with horizontal alignment, tilt, and second horizontal alignment
nk_deramp = xdem.coreg.NuthKaab() + xdem.coreg.Deramp() + xdem.coreg.NuthKaab()
nk_deramp.fit(dem_pleiades, reproj_dem_spot, inlier_mask=inlier_mask, verbose=True)
aligned_dem_spot = nk_deramp.apply(reproj_dem_spot)

# Save co-registered elevation differences to file
dh = dem_pleiades - aligned_dem_spot
dh.save('/home/atom/ongoing/work_stderr_dem/case_study_montblanc/dh_Pleiades-SPOT6_Mont-Blanc_NK_Deramp.tif')

# For Figure S2: saving independent steps of co-registration
nk = xdem.coreg.NuthKaab()
nk.fit(reproj_dem_spot, dem_pleiades, inlier_mask=inlier_mask, verbose=True)
aligned_dem_pleiades = nk.apply(dem_pleiades)

dh_nk = aligned_dem_pleiades - dem_pleiades
fn_nk = os.path.join(os.path.dirname(fn_dem_pleiades), 'dh_shift_nk_Pleiades.tif')
dh_nk.save(fn_nk)

deramp = xdem.coreg.Deramp()
deramp.fit(reproj_dem_spot, aligned_dem_pleiades, inlier_mask=inlier_mask, verbose=True)
deramped_dem_pleiades = deramp.apply(aligned_dem_pleiades)

dh_deramp = deramped_dem_pleiades - aligned_dem_pleiades
fn_deramp = os.path.join(os.path.dirname(fn_dem_pleiades), 'dh_shift_deramp_Pleiades.tif')
dh_deramp.save(fn_deramp)
