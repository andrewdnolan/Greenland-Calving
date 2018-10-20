"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Original method. Not optimized. Use NDWI_calc.py
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#!/usr/bin/env python
import os
import osr
import sys
import gdal
import argparse
import subprocess
import numpy as np
from argparse import RawTextHelpFormatter

def espg_code(gdal_ds):

    proj = osr.SpatialReference(wkt=gdal_ds.GetProjection())
    espg = int(proj.GetAttrValue('AUTHORITY',1))

    return espg

def reproject(raster_fn,t_srs=None):
    '''
        -   Right have it defaulting to Polar Stereo and manually inputting the parameters. Could and should change to use the EPSG code to alloow for
            other projections.

        -   Ideally wouldn't need to do this on the command line but due to errors with GDAL_DATA variable this is the sutiable alternative for right now
    '''
    no_f_type = os.path.splitext(raster_fn)[0]
    PS_fn = no_f_type + 'PS.TIF'

    if os.path.exists(PS_fn) == False:
        print ('Raster will be reprojected with Polar Stereo (ESPG: 3413)')
        #warp_cmd = "gdalwarp -t_srs '+proj=stere +lat_ts=70 +lat_0=90 +lon_0=-45+y_0=0 +x_0=0 +k=1 +datum=WGS84 +units=m'" + ' ' + raster_fn + ' ' + PS_fn
        warp_cmd = 'gdalwarp -t_srs "EPSG:' + str(t_srs) + '"'+ ' ' + raster_fn + ' ' + PS_fn
        subprocess.call(warp_cmd, shell=True)
    return PS_fn

def read_raster(raster,t_srs=None):
    ds = gdal.Open(raster)
    espg = espg_code(ds)

    if t_srs == None:
        t_srs=None
    else:
        if t_srs != espg:
            raster = reproject(raster,t_srs)
            ds = gdal.Open(raster)
            espg = espg_code(ds)

    ds_band = ds.GetRasterBand(1)
    band = ds_band.ReadAsArray()
    band.astype(np.float64)

    return band, ds

def calc_ndwi(Green, NIR, QA):
    valid_pixels = (QA !=1) & (Green + NIR != 0)

    NDWI= np.empty_like(Green, np.float64)
    NDWI[:] = -9999

    Green = Green[valid_pixels].astype(np.float64)
    NIR = NIR[valid_pixels].astype(np.float64)
    NDWI[valid_pixels] = (Green - NIR) / (Green + NIR)

    return NDWI

def ndwi_mask(NDWI,pix_thresh=0.15):
    valid_pixels = (NDWI != -9999) & (NDWI >= pix_thresh)

    NDWI_mask= np.empty_like(NDWI, np.float64)
    NDWI_mask[:] = 0

    NDWI = NDWI[valid_pixels].astype(np.float64)

    NDWI_mask[valid_pixels] = 1
    return NDWI_mask

def write_geotiff(Mask_array,raster_ds,input_raster_fn,output_raster_fn = None):
    #figuring out the name of the new file
    if output_raster_fn != None:
        mask_fn = output_raster_fn
    else:
        parent_dir, fn = os.path.split(input_raster_fn)
        mask_fn = parent_dir + '/'+'_'.join(fn.split('_')[:-1])+'_NDWI_mask.TIF'

    cols,rows = Mask_array.shape
    trans = raster_ds.GetGeoTransform()
    proj = raster_ds.GetProjection()
    outdriver = gdal.GetDriverByName("GTiff")
    outdata = outdriver.Create(mask_fn, rows, cols, 1, gdal.GDT_Float64)
    outdata.SetGeoTransform(trans)

    outdata.GetRasterBand(1).WriteArray(Mask_array)
    outdata.GetRasterBand(1).SetNoDataValue(0)
    outdata.SetProjection(proj)
    outdata.FlushCache()

    return mask_fn

#def replace_nan_w_0():


def getparser():
    description = ('Command line utility to calculate the Normailzed Differenced Water Index (NDWI) for a bulk order directory from earth Explorer. ')
    parser = argparse.ArgumentParser(description=description,formatter_class=RawTextHelpFormatter)
    parser.add_argument('src_dir',type=str,help='The source directory. Should be all of the available Landsat Scenes for a glacier. Should be a Bulk-Order directory from Earth Explorer')
    parser.add_argument('-t_srs',type=str, default = None, help='The target projection is the Landsat Scenes need to be reprojected')
    parser.add_argument('-pix_thres',type=float, default = 0.15, help = 'minimum value for the NDWI pixel value to be considered valid water.')
    parser.add_argument('-m',action='store_false',help='Whether the output should be a mask. Add flag to command line input to output mask.')
    return parser

def main():
    parser = getparser()
    args = parser.parse_args()

    src_dir = args.src_dir
    t_srs = args.t_srs
    pix_thres = args.pix_thres

    cwd = os.getcwd()

    for root,dirs,files in os.walk(src_dir):
        '''
        Needs to be sorted out. The inital walk through the directory with the files still zipped show up as file and not as dirs. But once they
        unzipped they will be show up as dirs which will complicate the os.walk process since they will be located in seperate list.
        '''

        # This is the inital walk for when the files are still zipped.
        for line, fn in enumerate(files):
            #compile a list of all of the scene directories and their file paths to itterate through
            if 'L' in fn:
                if fn.endswith('.tar.gz') == True and os.path.isdir(root + '/'+ fn.split('.')[0]) == False:
                    os.mkdir(root + '/' + fn.split('.')[0])
                    unzip_cmd = 'tar -xvzf ' + root+ '/' + fn + ' -C ' + root + '/' + fn.split('.')[0]
                    subprocess.call(unzip_cmd,shell=True)


        # This is the second walk for once the files are unzipped and appear as dirs.
        for line, dn in enumerate(dirs):
            #Might need to add more conditionals here based on what the naming convention is for all of Landsat 8
            if dn.startswith('LC') == True:
                bands = os.listdir(root + '/' +dn)
                for i,band in enumerate(bands):
                    if band.endswith('B3.TIF') == True:
                        green_fn = cwd + '/' + root + '/' + dn + '/' + bands[i]
                        green_raster, green_ds = read_raster(green_fn, t_srs=t_srs)
                    elif band.endswith('B5.TIF') == True:
                        nir_fn = cwd + '/' + root + '/' + dn + '/' + bands[i]
                        nir_raster, nir_ds = read_raster(nir_fn, t_srs=t_srs)
                    elif band.endswith('BQA.TIF') == True:
                        qa_fn = cwd + '/' + root + '/' + dn + '/' + bands[i]
                        qa_raster, qa_ds = read_raster(qa_fn, t_srs=t_srs)

                NDWI = calc_ndwi(green_raster,nir_raster,qa_raster)

                # Need to add the fucntionality to choose whetehr the ooutput is a mask or just the NDWI.
                NDWI_mask = ndwi_mask(NDWI,pix_thres)

                mask_fn = write_geotiff(NDWI_mask,green_ds,green_fn)



if __name__ == '__main__':
    main()
