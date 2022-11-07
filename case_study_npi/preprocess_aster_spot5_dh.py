"""Pre-process ASTER and SPOT-5 DEMs into elevation changes for the NPI case study"""
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    import xdem
import geoutils as gu
import numpy as np
import matplotlib.pyplot as plt

# Open DEMs and outlines
fn_glacier_outlines = '/home/atom/data/inventory_products/RGI/00_rgi60_neighb_merged/17_rgi60_SouthernAndes/17_rgi60_SouthernAndes.shp'
fn_dem_spot = '/home/atom/ongoing/work_stderr_dem/case_study_npi/SPOT5_2012-03-18_NPI_NDV.tif'
fn_dem_aster = '/home/atom/ongoing/work_stderr_dem/case_study_npi/AST_L1A_00303182012144228/ASTER_NPI_00303182012144228_Z.tif'
fn_corr='/home/atom/ongoing/work_stderr_dem/case_study_npi/AST_L1A_00303182012144228/ASTER_NPI_00303182012144228_CORR.tif'

dem_spot = gu.Raster(fn_dem_spot)
dem_aster = xdem.DEM(fn_dem_aster)
glacier_outlines = gu.Vector(fn_glacier_outlines)

# Reproject on the ASTER DEM
reproj_dem_spot = dem_spot.reproject(dem_aster)
init_dh = reproj_dem_spot - dem_aster

# Open quality of stereo-correlation raster, and remove low corrections for alignment
corr = gu.Raster(fn_corr)
mask_highcorr = corr.data>70.
mask_noglacier = ~glacier_outlines.create_mask(dem_aster)
mask_nooutliers = np.abs(init_dh.data - np.nanmedian(init_dh.data))<3*xdem.spatialstats.nmad(init_dh.data)

# Create inlier mask
inlier_mask = np.logical_and.reduce((mask_noglacier, mask_nooutliers, mask_highcorr))

# Coregistration pipeline with horizontal alignment, tilt, and second horizontal alignment
nk_deramp = xdem.coreg.NuthKaab() + xdem.coreg.Deramp() + xdem.coreg.NuthKaab()
nk_deramp.fit(dem_aster, reproj_dem_spot, inlier_mask=inlier_mask, verbose=True)
aligned_dem_spot = nk_deramp.apply(reproj_dem_spot)

# Save co-registered elevation differences to file
dh = dem_aster - aligned_dem_spot
dh.save('/home/atom/ongoing/work_stderr_dem/case_study_npi/dh_ASTER-SPOT5_NPI_NK_Deramp.tif')