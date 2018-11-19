#!/usr/bin/env python

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
v2. is a rasterio based for the geospatial process, vs the GDAL Based v1.


Command line utility for calculating the NDWI for a input directory of bulk-Order
Landsat scenes from Earth Explorerself.

Another option may be to parrallelize the decompression of the .tar.gz Landsat
scenes. This will present a significant slow down if it is only computed by one
core. This will especially be the case when the script is executed on an direcotry
containg 50-100 Landsat scenes to be processed.
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
import os
import sys
import argparse
import rasterio
import subprocess
import numpy as np
import numpy.ma as ma
from rio_toa import toa_utils
from rio_toa import reflectance
from argparse import RawTextHelpFormatter

#Get a list of the unzipped Landat directories to process
def walk_dir(src_dir):
    dir_list = []
    #Checks if the zipped file has been unzipped. If not, unzips
    for root,dirs,files in os.walk(src_dir):
        for line, fn in enumerate(files):
            if 'L' in fn:
                if fn.endswith('.tar.gz') == True and os.path.isdir(root + '/'+ fn.split('.')[0]) == False:
                    os.mkdir(root + '/' + fn.split('.')[0])
                    unzip_cmd = 'tar -xvzf ' + root+ '/' + fn + ' -C ' + root + '/' + fn.split('.')[0]
                    subprocess.call(unzip_cmd,shell=True)
    # Second walk through after the unzipping (if neccessary)
    for root,dirs,files in os.walk(src_dir):
        for line, dn in enumerate(dirs):
            #Might need to add more conditionals here based on what the naming convention is for all of Landsat 8
            if dn.startswith('LC') == True or dn.startswith('LT') == True:
                Ldir = root + '/' + dn
                dir_list.append(Ldir)
    return list(set(dir_list))

# Get source EPSG code to test if input secenes need to reporjected
def epsg_code(raster_fn):
    with rasterio.open(raster_fn, 'r') as src:
        epsg = src.crs.to_epsg()
    return epsg

#TOA reflectace corection for Landsat 8
def TOA_corc(raster_fn):
    processes = 4
    dtype = 'uint16'
    creation_options = {'nodata':0, 'compress':'deflate', 'predict':2}
    rescale_factor = 55000

    scene_mtl = "/".join(raster_fn.split('/')[0:-1]) +'/'+ '_'.join(raster_fn.split('/')[-1].split('_')[0:7]) + '_MTL.txt'
    band_num = raster_fn.split('/')[-1].split('_')[7][1]
    if os.path.exists(scene_mtl) == False:
        sys.exit('Can not find Meta Data file (MTL.txt) for scene')

    output_fp = "/".join(raster_fn.split('/')[0:-1]) + '/' + '_'.join(raster_fn.split('/')[-1].split('_')[0:7]) + '_TOA_' + '_'.join(raster_fn.split('/')[-1].split('_')[7:])
    if not os.path.exists(output_fp) and not os.path.exists(os.path.splitext(output_fp)[0] + '_30PS.TIF'):
        reflectance.calculate_landsat_reflectance([raster_fn], scene_mtl, output_fp, rescale_factor=rescale_factor, creation_options=creation_options, bands=[band_num], dst_dtype=dtype, processes=processes, pixel_sunangle=True)

#reproject and resample the input Landsat Scenes
def reproject(raster_fn,t_srs=3413):
    '''
        -   This method changes in the input resolution in order to make them
            align with GIMP masks correctly

        -   Right have it defaulting to Polar Stereo and manually inputting the
            parameters. Could and should change to use the EPSG code to allow for
            other projections.

    '''
    no_f_type = os.path.splitext(raster_fn)[0]
    PS_fn = no_f_type + '_30PS.TIF'

    if os.path.exists(PS_fn) == False:
        print ('Raster will be reprojected with Polar Stereo (EPSG: 3413)')
        warp_cmd = 'gdalwarp -t_srs "EPSG:' + str(t_srs) + '"'+ ' ' + raster_fn + ' ' + PS_fn + ' -tr 30 30'
        subprocess.call(warp_cmd, shell=True)
    return PS_fn

#determine the intersection of the Mask and Landsat scene
def interection(scene_1, scene_2):
    '''
    rasterio based function to determine the overlap between two rasters
    Returns a tuple of the bounds of intersection which can be used by
    rasterio.window functionality
    '''
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

#Read in the rasters (B3, B5, GIMP_mask) with rasterio
def read_raster(land_dir, GIMP_fn):
    bands = os.listdir(land_dir)
    for i,band in enumerate(bands):
        if band.endswith('B3_30PS.TIF') == True:
            g_fn = land_dir + '/' + bands[i]
        elif band.endswith('B5_30PS.TIF') == True:
            nir_fn = land_dir + '/' + bands[i]
        elif band.endswith('B6_30PS.TIF') == True:
            swnir_fn = land_dir + '/' + bands[i]

    interesction = interection(g_fn, GIMP_fn)
    '''
    Should set up a break in the script here if the landsat scene and the
    GIMP mask don't overlap
    '''
    with rasterio.open(g_fn) as src:
        sub_window = src.window(*interesction, precision = 15)
        g_bounds = src.bounds
        g = src.read(1, window = sub_window)

    with rasterio.open(nir_fn) as src:
        sub_window = src.window(*interesction, precision = 15)
        nir= src.read(1, window = (sub_window))

    # with rasterio.open(qa_fn) as src:
    #     sub_window = src.window(*hm_interesction)
    #     qa = src.read(1, window = sub_window)

    with rasterio.open(GIMP_fn) as src:
        sub_window = src.window(*interesction, precision = 15)
        GIMP_bounds = src.bounds
        gimp_trans = src.transform
        GIMP_mask = src.read(1, window = sub_window)

    src_att = [src.crs.to_epsg(), sub_window, src.transform]
    return g, nir, GIMP_mask, src_att

#This ain't pretty but has to get done
def inner_bounds (LS_array, GIMP_array):
    if GIMP_array.shape[0] > LS_array.shape[0]:
        inx = LS_array.shape[0]
    elif LS_array.shape[0] > GIMP_array.shape[0]:
        inx = GIMP_array.shape[0]
    elif LS_array.shape[0] == GIMP_array.shape[0]:
        #this doesn't matter, just need to assign the var
        inx = LS_array.shape[0]

    if GIMP_array.shape[1] > LS_array.shape[1]:
        iny = LS_array.shape[1]
    elif LS_array.shape[1] > GIMP_array.shape[1]:
        iny = GIMP_array.shape[1]
    elif LS_array.shape[1] == GIMP_array.shape[1]:
        #this doesn't matter, just need to assign the var
        iny = LS_array.shape[1]

    return inx, iny

#NDWI calculation
def NDWI_calc(Green, NIR):
    NDWI = (Green.astype(float) - NIR.astype(float)) / (Green + NIR)

    return NDWI

#Mask generator from NDWI array based on given pixel threshold
def NDWI_npmask (NDWI_array, pix_thres):
    NDWI_array = np.ma.masked_where( NDWI_array < pix_thres, NDWI_array)
    NDWI_mask = NDWI_array.filled(fill_value = 0)
    NDWI_mask = NDWI_mask > 0
    NDWI_mask = NDWI_mask.astype(np.uint8)

    return (NDWI_mask)

def write_mask(array,land_dir, pix_thes, src_att):
    template = [fn for fn in os.listdir(land_dir) if fn.endswith('MTL.txt')][0]
    mask_fn = land_dir +'/' + '_'.join(template.split('_')[0:-1]) + '_{}MSK.TIF'.format(pix_thes)

    with rasterio.open(mask_fn, 'w', driver = 'GTiff', height = array.shape[0], width = array.shape[1], count = 1, dtype = rasterio.uint8, crs='EPSG:{}'.format(src_att[0]), transform = rasterio.windows.transform(src_att[1], src_att[2])) as dst:
        dst.write(array, 1)

def getparser():
    description = ('Command line utility to calculate the Normailzed Differenced Water Index (NDWI) for a bulk order directory from earth Explorer. ')
    parser = argparse.ArgumentParser(description=description,formatter_class=RawTextHelpFormatter)
    parser.add_argument('src_dir',type=str,help='The source directory. Should be all of the available Landsat Scenes for a glacier. Should be a Bulk-Order directory from Earth Explorer (Ex. "Bulk_Order_930842")')
    parser.add_argument('pix_thres',type=float, default = 0.15, help = 'minimum value for the NDWI pixel value to be considered valid water.')
    parser.add_argument('GIMP_fn', default = None, help = 'GIMP tile filename for seclected Glacier')
    parser.add_argument('-t_srs',type=str, default = 3413, help='The target projection is the Landsat Scenes need to be reprojected')

    return parser

def main():
    parser = getparser()
    args = parser.parse_args()

    main_dir = args.src_dir
    t_srs = args.t_srs
    pix_thres = args.pix_thres
    GIMP_fn = args.GIMP_fn

    Land_dirs = walk_dir(main_dir)

    #check if the mask is projected right
    if t_srs != epsg_code(GIMP_fn):
        GIMP_fn = reproject(GIMP_fn)

    #this is incase you want to do the claculations on just one Landsat scene
    if Land_dirs == []:
        for band in os.listdir(main_dir):
            if band.endswith('T1_B3.TIF') or band.endswith('T1_B5.TIF'):
                TOA_corc(main_dir + '/' + band)
        for band in os.listdir(main_dir):
            if band.endswith('TOA_B3.TIF') or band.endswith('TOA_B5.TIF'):
                if t_srs != epsg_code(main_dir + '/' + band):
                    reproject(main_dir + '/' + band)

        green, nir, GIMP, src_att = read_raster(main_dir, GIMP_fn)

        row_col = inner_bounds(green, GIMP)

        g_m = ma.array(green[:row_col[0],:row_col[1]], mask = (GIMP == 0), dtype =np.float64)
        g_m = ma.array(g_m, mask = (1 >= g_m), dtype =np.float64)

        nir_m = ma.array(nir[:row_col[0],:row_col[1]], mask = (GIMP == 0), dtype =np.float64)
        nir_m = ma.array(nir_m, mask = (1 >= nir_m), dtype =np.float64)

        ndwi_array = NDWI_calc(g_m, nir_m)

        ndwi_m_array = NDWI_npmask(ndwi_array, 0.25)

        write_mask(ndwi_m_array, main_dir, 0.25, src_att)

    if Land_dirs != []:
        for src_dir in Land_dirs:
            for band in os.listdir(src_dir):
                if band.endswith('T1_B3.TIF') or band.endswith('T1_B5.TIF') or band.endswith('T2_B5.TIF') or band.endswith('T2_B3.TIF'):
                    TOA_corc(src_dir + '/' + band)
            for band in os.listdir(src_dir):
                if band.endswith('TOA_B3.TIF') or band.endswith('TOA_B5.TIF'):
                    if t_srs != epsg_code(src_dir + '/' + band):
                        reproject(src_dir + '/' + band)
            green, nir, GIMP, src_att = read_raster(src_dir, GIMP_fn)

            row_col = inner_bounds(green, GIMP)

            g_m = ma.array(green[:row_col[0],:row_col[1]], mask = (GIMP == 0), dtype =np.float64)
            g_m = ma.array(g_m, mask = (1 >= g_m), dtype =np.float64)

            nir_m = ma.array(nir[:row_col[0],:row_col[1]], mask = (GIMP == 0), dtype =np.float64)
            nir_m = ma.array(nir_m, mask = (1 >= nir_m), dtype =np.float64)

            ndwi_array = NDWI_calc(g_m, nir_m)

            ndwi_m_array = NDWI_npmask(ndwi_array, 0.25)

            write_mask(ndwi_m_array, src_dir, 0.25, src_att)

if __name__ == '__main__':
    main()
