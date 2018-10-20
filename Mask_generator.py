import os
#import osr
import sys
import rasterio
import numpy as np
import numpy.ma as ma

def interection(scene_1, scene_2):
    with rasterio.open(scene_1) as src:
        S1_bounds = src.bounds
    with rasterio.open(scene_2) as src:
        S2_bounds = src.bounds

    inter_tuple = []
    for i in range(len(tuple(S1_bounds))):
        if i == 0:
            if S1_bounds[i] < S2_bounds[i]:
                inter_tuple.append(S2_bounds[i])
            else:
                inter_tuple.append(S1_bounds[i])
        if i ==1:
            if S1_bounds[i] < S2_bounds[i]:
                inter_tuple.append(S2_bounds[i])
            else:
                inter_tuple.append(S1_bounds[i])
        if i == 2:
            if S1_bounds[i] > S2_bounds[i]:
                inter_tuple.append(S2_bounds[i])
            else:
                inter_tuple.append(S1_bounds[i])
        if i == 3:
            if S1_bounds[i] > S2_bounds[i]:
                inter_tuple.append(S2_bounds[i])
            else:
                inter_tuple.append(S1_bounds[i])
    inter_tuple = tuple(inter_tuple)

    return(inter_tuple)

# src_dir = '/Users/glaciologygroup/Greenland_Calving/Landsat_NDWI/Sample_data/Bulk_Order_931485/Landsat_8_OLI_TIRS_C1_Level-1'
# #mask_fn =
# mask_dir = '/Users/glaciologygroup/Greenland_Calving/Landsat_NDWI/GIMP_masks/Landsat_Masks'
# for sroot,sdirs,sfiles in os.walk(src_dir):
#     '''
#     Needs to be sorted out. The inital walk through the directory with the files still zipped show up as file and not as dirs. But once they
#     unzipped they will be show up as dirs which will complicate the os.walk process since they will be located in seperate list.
#     '''
#
#     # This is the inital walk for when the files are still zipped.
#     # Could be a  good idea to parallelize this. Multithreading this and giving each core a direcotry to unzip might not be the worst idea
#     # especially when you will have 100+ scenes per direcotry of each glacier.
#     count = 0
#     for line, fn in enumerate(sfiles):
#         #compile a list of all of the scene directories and their file paths to itterate through
#         if 'L' in fn:
#             #Test if all the scenes are unzipped. Will need them to be unzipped
#             #since we will need the extent of the raster to clip the mask too
#             if fn.endswith('.tar.gz') == True and os.path.isdir(sroot + '/'+ fn.split('.')[0]) == False:
#                 x =1
#                 #os.mkdir(root + '/' + fn.split('.')[0])
#                 #unzip_cmd = 'tar -xvzf ' + root+ '/' + fn + ' -C ' + root + '/' + fn.split('.')[0]
#                 #subprocess.call(unzip_cmd,shell=True)
#
# pr_list = []
# fn_list = []
# dic = {}
# #now rewalk direcotry to get the new files if the direcotry was just unzipped
# for sroot,sdirs,sfiles in os.walk(src_dir):
#     for line, dn in enumerate(sdirs):
#         if dn.startswith('LC'):
#             path_row = dn.split('_')[2]
#             pr_list.append(path_row)
#
# pr_list = list(set(pr_list))
# for i in pr_list:
#     dic[str(i)] = []
#
# for sroot,sdirs,sfiles in os.walk(src_dir):
#     for line, fn in enumerate(sfiles):
#         for key in dic:
#             if key in fn and fn.endswith('PS.TIF') and not dic[key]:
#                 dic[key].append(sroot + '/' +fn)
# for key in dic:
#     data = gdal.Open(dic[key][0], GA_ReadOnly)
#     geoTransform = data.GetGeoTransform()
#     minx = geoTransform[0]
#     maxy = geoTransform[3]
#     maxx = minx + geoTransform[1] * data.RasterXSize
#     miny = maxy + geoTransform[5] * data.RasterYSize
#
# rpjected_fn = '/Users/glaciologygroup/Greenland_Calving/Landsat_NDWI/Sample_data/Bulk_Order_931485/Landsat_8_OLI_TIRS_C1_Level-1/LC08_L1TP_009011_20180730_20180730_01_RT/LC08_L1TP_009011_20180730_20180730_01_RT_B5PS.TIF'
#
# mask_dir = '/Users/glaciologygroup/Greenland_Calving/Landsat_NDWI/GIMP_masks/Landsat_Masks/'
# rpjected_fn = '/Users/glaciologygroup/Greenland_Calving/Landsat_NDWI/Harold_Moltke/Bulk_Order_948864/Landsat_8_OLI_TIRS_C1_Level-1/LC08_L1TP_031005_20150630_20170407_01_T1/LC08_L1TP_031005_20150630_20170407_01_T1_B5_30PS.TIF'
# tile_fn = '/Users/glaciologygroup/Greenland_Calving/Landsat_NDWI/GIMP_masks/GimpIceMask_30m_tile0_4_v1.1_PS.tif'
#
# def mask_existance(rpjected_fn, mask_dir):
#     mask_fns = [os.path.join(path, filename) for path, dirnames, filenames in os.walk(mask_dir) for filename in filenames if filename.endswith('.tif')]
#     path_row = rpjected_fn.split('/')[-1].split('_')[2]
#
#     for fn in mask_fns:
#         if path_row in fn:
#             exists = True
#         else:
#             exists = False
#
#     if not mask_fns: #for when the mask dir is empty
#         exists = False
#     return exists
#
# def mask_create(rpjected_fn, tile_fn, mask_dir):
#     path_row = rpjected_fn.split('/')[-1].split('_')[2]
#     mask_fn = mask_dir + 'GIMP_LS_' + path_row + '.tif'
#     pr_shp = mask_dir + 'GIMP_LS_' + path_row + '.shp'
#     data = gdal.Open(rpjected_fn, GA_ReadOnly)
#     geoTransform = data.GetGeoTransform()
#     minx = geoTransform[0]
#     maxy = geoTransform[3]
#     maxx = minx + geoTransform[1] * data.RasterXSize
#     miny = maxy + geoTransform[5] * data.RasterYSize
#     clip_cmd = 'gdal_translate -projwin ' + ' '.join([str(x) for x in [minx, maxy, maxx, miny]]) + ' -of GTiff ' + rpjected_fn + ' ' + mask_fn
#     subprocess.call(clip_cmd,shell=True)
#     return mask_fn
#
# mask_check_op = mask_check(rpjected_fn, mask_dir)
#
# if mask_check_op == False:
#     test_mask = mask_create(rpjected_fn, tile_fn, mask_dir)
#     print(test_mask)
# #landsat_fn ='/Users/glaciologygroup/Greenland_Calving/Landsat_NDWI/Harold_Moltke/Bulk_Order_948864/Landsat_8_OLI_TIRS_C1_Level-1/LC08_L1TP_029006_20130930_20170502_01_T1/LC08_L1TP_029006_20130930_20170502_01_T1_B5_PS.TIF'
# Binary_Mask = '/Users/glaciologygroup/Greenland_Calving/Landsat_NDWI/GIMP_masks/Landsat_Masks/GIMP_LS_031005.tif'
#
# with rasterio.open(rpjected_fn) as src:
#     src_crs = src.crs
#     src_transform = src.transform
#     src_bounds = src.bounds
#     src_width = src.width
#     src_height = src.height
#     Landsat_array = src.read(1)
#
#
#
#
# with rasterio.open(Binary_Mask) as src:
#     mask = src.read(1)
#     mask_bounds = src.bounds
#
# Masked_Landsat = np.ma.array(Landsat_array, mask = mask!=0)
