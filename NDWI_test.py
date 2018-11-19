#########################################################################################################
#this script is to trouble shoot problems with the NDWI calcuitons
# PROBLEMS NEEDING TROUBLE SHOOTING:
#       - make it a binary output of 0 &1. Currently can only have successful outputs if the no data value is
#         set to 0, but then the raster is filled with nans and 1's
#       - Make this faster!!!
#           - we have tried reading in the raster with blocks which is supposed to be more effeicnet
#           - could try to parallelize it through mutli proccess
#               - if i try to multi thread it proabably will run in to GIL errors pretty quickly
##########################################################################################################
import os
os.chdir('/Users/glaciologygroup/Greenland_Calving/Landsat_NDWI/NDWI_Func')
from NDWI_calc_v2 import *

src_dir = '/Users/glaciologygroup/Greenland_Calving/Landsat_NDWI/thresh_test/LC08_L1TP_017008_20140828_20170420_01_T1'
gimp_mask = '/Users/glaciologygroup/Greenland_Calving/Landsat_NDWI/GIMP_masks/GimpIceMask_30m_tile1_3_v1.1_PS.tif'


for band in os.listdir(src_dir):
    if band.endswith('B3.TIF') or band.endswith('B5.TIF'):
        TOA_corc(src_dir + '/' + band)
for band in os.listdir(src_dir):
    if band.endswith('TOA_B3.TIF') or band.endswith('TOA_B5.TIF'):
        reproject(src_dir + '/' + band)

green, nir, GIMP, src_att = read_raster(src_dir, gimp_mask)

row_col = inner_bounds(green, GIMP)

g_m = ma.array(green[:row_col[0],:row_col[1]], mask = (GIMP == 0), dtype =np.float64)
g_m = ma.array(g_m, mask = (1 >= g_m), dtype =np.float64)

nir_m = ma.array(nir[:row_col[0],:row_col[1]], mask = (GIMP == 0), dtype =np.float64)
nir_m = ma.array(nir_m, mask = (1 >= nir_m), dtype =np.float64)

ndwi_array = NDWI_calc(g_m, nir_m)

ndwi_m_array = NDWI_npmask(ndwi_array, 0.25)

write_mask(ndwi_m_array, src_dir, 0.25, src_att)
