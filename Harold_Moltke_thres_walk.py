import os
import sys
import rasterio
import numpy as np
import pandas as pd
import numpy.ma as ma
import matplotlib.pyplot as plt
from NDWI_calc import calc_ndwi, ndwi_mask
from Mask_generator import interection

def NDWI_calc(Green, NIR):
    NDWI = (Green.astype(float) - NIR.astype(float)) / (Green + NIR)

    return NDWI

def NDWI_npmask (NDWI_array, pix_thres):
    NDWI_array = np.ma.masked_where( NDWI_array < pix_thres, NDWI_array)
    NDWI_mask = NDWI_array.filled(fill_value = 0)
    NDWI_mask = NDWI_mask > 0
    NDWI_mask = NDWI_mask.astype(np.uint8)

    return (NDWI_mask)

#http://karthur.org/2015/clipping-rasters-in-python.html

# ######################### Walk through for Harold Moltke ######################
# hm_dir = '/Users/glaciologygroup/Greenland_Calving/Landsat_NDWI/Harold_Moltke/Bulk_Order_948864/Landsat_8_OLI_TIRS_C1_Level-1/LC08_L1TP_031005_20130827_20170502_01_T1'
# HM_GIMP = '/Users/glaciologygroup/Greenland_Calving/Landsat_NDWI/GIMP_masks/GimpIceMask_30m_tile0_4_v1.1_PS.tif'
#
# bands = os.listdir(hm_dir)
# for i,band in enumerate(bands):
#     if band.endswith('B3_30PS.TIF') == True:
#         hm_g_fn = hm_dir + '/' + bands[i]
#     elif band.endswith('B5_30PS.TIF') == True:
#         hm_nit_fn = hm_dir + '/' + bands[i]
#     # elif band.endswith('BQA_30PS.TIF') == True:
#     #     qa_fn = hm_dir + '/' + bands[i]
#
# hm_interesction = interection(hm_g_fn, HM_GIMP)
#
# with rasterio.open(hm_g_fn) as src:
#     sub_window = src.window(*hm_interesction, precision = 15)
#     hm_g_bounds = src.bounds
#     hm_g = src.read(1, window = sub_window)
#
# with rasterio.open(hm_nit_fn) as src:
#     sub_window = src.window(*hm_interesction, precision = 15)
#     hm_nir= src.read(1, window = (sub_window))
#
# # with rasterio.open(qa_fn) as src:
# #     sub_window = src.window(*hm_interesction)
# #     qa = src.read(1, window = sub_window)
#
# with rasterio.open(HM_GIMP) as src:
#     sub_window = src.window(*hm_interesction, precision = 15)
#     hm_GIMP_bounds = src.bounds
#     hm_gimp_trans = src.transform
#     GIMP_mask = src.read(1, window = sub_window)
#
# # if hm_g_bounds != hm_GIMP_bounds:
# #     sys.exit('Need to check source rasters. Mask and Landsat Scene do not align Properly')
#
# hm_g_m = ma.array(hm_g[:GIMP_mask.shape[0],:GIMP_mask.shape[1]], mask = (GIMP_mask == 0), dtype =np.float64)
# hm_g_m = ma.array(hm_g_m, mask = (1 >= hm_g_m), dtype =np.float64)
# #hm_g= hm_g_m.filled(fill_value=-9999)
#
# hm_nir_m = ma.array(hm_nir[:GIMP_mask.shape[0],:GIMP_mask.shape[1]], mask = (GIMP_mask == 0), dtype =np.float64)
# hm_nir_m = ma.array(hm_nir_m, mask = (1 >= hm_nir_m), dtype =np.float64)
# #hm_nir= hm_nir_m.filled(fill_value=-9999)
#
# thresh_parm_list = [i / 100.0 for i in range(0, 105, 5)]
#
# hm_ndwi = NDWI_calc(hm_g_m, hm_nir_m)
# hm_ndwi_masks = []
# hm_mask_areas = []
#
# for i in thresh_parm_list:
#     thresh_array = NDWI_npmask(hm_ndwi, i)
#     thresh_unique = np.unique(thresh_array, return_counts=True)
#     hm_ndwi_masks.append(thresh_array)
#     if len(thresh_unique[0]) != 2:
#         if thresh_unique[0][0] == 0:
#             area = 0
#             hm_mask_areas.append(area)
#     else:
#         area = (np.unique(thresh_array, return_counts=True)[1][1] * 30.0 * 30.0) / (1000.0**2)
#         hm_mask_areas.append(area)
#         mask_fn = '{}/HM_{}_thresh.tif'.format(hm_dir, i)
#         with rasterio.open(mask_fn, 'w', driver = 'GTiff', height = thresh_array.shape[0], width = thresh_array.shape[1], count = 1, dtype = rasterio.uint8, crs='EPSG:3413', transform = rasterio.windows.transform(sub_window, src.transform)) as dst:
#             dst.write(thresh_array, 1)
#
#  
#
# area_df = pd.DataFrame({'HM':hm_mask_areas}, index = thresh_parm_list)


# ################ Upernavik N / S walk through ################################
# up_dir = '/Users/glaciologygroup/Greenland_Calving/Landsat_NDWI/Upernavik_N_S/Bulk_Order_931150/Landsat_8_OLI_TIRS_C1_Level-1/LC08_L1TP_017008_20140828_20170420_01_T1'
# up_GIMP = '/Users/glaciologygroup/Greenland_Calving/Landsat_NDWI/GIMP_masks/GimpIceMask_30m_tile1_3_v1.1_PS.tif'
#
# bands = os.listdir(up_dir)
# for i,band in enumerate(bands):
#     if band.endswith('B3_30PS.TIF') == True:
#         up_g_fn = up_dir + '/' + bands[i]
#     elif band.endswith('B5_30PS.TIF') == True:
#         up_nit_fn = up_dir + '/' + bands[i]
#     # elif band.endswith('BQA_30PS.TIF') == True:
#     #     qa_fn = up_dir + '/' + bands[i]
#
# up_interesction = interection(up_g_fn, up_GIMP)
#
# with rasterio.open(up_g_fn) as src:
#     sub_window = src.window(*up_interesction)
#     up_g_bounds = src.bounds
#     up_g = src.read(1, window = sub_window)
#
# with rasterio.open(up_nit_fn) as src:
#     sub_window = src.window(*up_interesction)
#     up_nir= src.read(1, window = sub_window)
#
# # with rasterio.open(qa_fn) as src:
# #     sub_window = src.window(*up_interesction)
# #     qa = src.read(1, window = sub_window)
#
# with rasterio.open(up_GIMP) as src:
#     sub_window = src.window(*up_interesction)
#     up_GIMP_bounds = src.bounds
#     GIMP_mask = src.read(1, window = sub_window)
#
# # if up_g_bounds != up_GIMP_bounds:
# #     sys.exit('Need to check source rasters. Mask and Landsat Scene do not align Properly')
#
# up_g_m = ma.array(up_g[:GIMP_mask.shape[0],:GIMP_mask.shape[1]], mask = (GIMP_mask == 0), dtype =np.float64)
# up_g_m = ma.array(up_g_m, mask = (1 >= up_g_m), dtype =np.float64)
# #up_g= up_g_m.filled(fill_value=-9999)
#
# up_nir_m = ma.array(up_nir[:GIMP_mask.shape[0],:GIMP_mask.shape[1]], mask = (GIMP_mask == 0), dtype =np.float64)
# up_nir_m = ma.array(up_nir_m, mask = (1 >= up_nir_m), dtype =np.float64)
# #up_nir= up_nir_m.filled(fill_value=-9999)
#
# up_ndwi = NDWI_calc(up_g_m, up_nir_m)
# up_ndwi_masks = []
# up_mask_areas = []
#
# for i in thresh_parm_list:
#     thresh_array = NDWI_npmask(up_ndwi, i)
#     thresh_unique = np.unique(thresh_array, return_counts=True)
#     up_ndwi_masks.append(thresh_array)
#     if len(thresh_unique[0]) != 2:
#         if thresh_unique[0][0] == 0:
#             area = 0
#             up_mask_areas.append(area)
#     else:
#         area = (np.unique(thresh_array, return_counts=True)[1][1] * 30.0 * 30.0) / (1000.0**2)
#         up_mask_areas.append(area)
#
#     print('Pix threshold - {}  Area - {}'.format(i, area))
#
# area_df['up'] = up_mask_areas
#
# ################### Walk through for Heimdal #################################
# hd_dir = '/Users/glaciologygroup/Greenland_Calving/Landsat_NDWI/Heimdal/LC08_L1TP_233016_20170712_20170726_01_T1'
# hd_GIMP = '/Users/glaciologygroup/Greenland_Calving/Landsat_NDWI/GIMP_masks/GimpIceMask_30m_tile2_0_and_3_0_v1.1_PS.tif'
#
# bands = os.listdir(hd_dir)
# for i,band in enumerate(bands):
#     if band.endswith('B3_30PS.TIF') == True:
#         hd_g_fn = hd_dir + '/' + bands[i]
#     elif band.endswith('B5_30PS.TIF') == True:
#         hd_nit_fn = hd_dir + '/' + bands[i]
#     # elif band.endswith('BQA_30PS.TIF') == True:
#     #     qa_fn = hd_dir + '/' + bands[i]
#
# hd_interesction = interection(hd_g_fn, hd_GIMP)
#
# with rasterio.open(hd_g_fn) as src:
#     sub_window = src.window(*hd_interesction)
#     hd_g_bounds = src.bounds
#     hd_g = src.read(1, window = sub_window)
#
# with rasterio.open(hd_nit_fn) as src:
#     sub_window = src.window(*hd_interesction)
#     hd_nir= src.read(1, window = sub_window)
#
# # with rasterio.open(qa_fn) as src:
# #     sub_window = src.window(*hd_interesction)
# #     qa = src.read(1, window = sub_window)
#
# with rasterio.open(hd_GIMP) as src:
#     sub_window = src.window(*hd_interesction)
#     hd_GIMP_bounds = src.bounds
#     GIMP_mask = src.read(1, window = sub_window)
#
#
# hd_g_m = ma.array(hd_g[:GIMP_mask.shape[0],:GIMP_mask.shape[1]], mask = (GIMP_mask == 0), dtype =np.float64)
# hd_g_m = ma.array(hd_g_m, mask = (1 >= hd_g_m), dtype =np.float64)
# #hd_g= hd_g_m.filled(fill_value=-9999)
#
# hd_nir_m = ma.array(hd_nir[:GIMP_mask.shape[0],:GIMP_mask.shape[1]], mask = (GIMP_mask == 0), dtype =np.float64)
# hd_nir_m = ma.array(hd_nir_m, mask = (1 >= hd_nir_m), dtype =np.float64)
# #hd_nir= hd_nir_m.filled(fill_value=-9999)
#
# hd_ndwi = NDWI_calc(hd_g_m, hd_nir_m)
# hd_ndwi_masks = []
# hd_mask_areas = []
#
# for i in thresh_parm_list:
#     thresh_array = NDWI_npmask(hd_ndwi, i)
#     thresh_unique = np.unique(thresh_array, return_counts=True)
#     hd_ndwi_masks.append(thresh_array)
#     if len(thresh_unique[0]) != 2:
#         if thresh_unique[0][0] == 0:
#             area = 0
#             hd_mask_areas.append(area)
#     else:
#         area = (np.unique(thresh_array, return_counts=True)[1][1] * 30.0 * 30.0) / (1000.0**2)
#         hd_mask_areas.append(area)
#
#     print('Pix threshold - {}  Area - {}'.format(i, area))
#
# area_df['hd'] = hd_mask_areas
# area_csv_fn = '/Users/glaciologygroup/Greenland_Calving/Landsat_NDWI/thresh_areas_csv.text'
# area_df.to_csv(area_csv_fn, sep=';')
