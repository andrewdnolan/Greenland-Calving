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
import osr
import gdal
import subprocess
import numpy as np
from NDWI_Func_block_io import *
import matplotlib.pyplot as plt

green_fn  = '/Users/glaciologygroup/Desktop/NDWI_Test/Bulk_Order_931485/Landsat_8_OLI_TIRS_C1_Level-1/LC08_L1TP_009011_20180730_20180730_01_RT/LC08_L1TP_009011_20180730_20180730_01_RT_B3.TIF'
nir_fn = '/Users/glaciologygroup/Desktop/NDWI_Test/Bulk_Order_931485/Landsat_8_OLI_TIRS_C1_Level-1/LC08_L1TP_009011_20180730_20180730_01_RT/LC08_L1TP_009011_20180730_20180730_01_RT_B5.TIF'
qa_fn = '/Users/glaciologygroup/Desktop/NDWI_Test/Bulk_Order_931485/Landsat_8_OLI_TIRS_C1_Level-1/LC08_L1TP_009011_20180730_20180730_01_RT/LC08_L1TP_009011_20180730_20180730_01_RT_BQA.TIF'

whole_array = read_n_write(green_fn, nir_fn, qa_fn, t_srs=3413)
plt.imshow(whole_array, clim=(0,1), cmap = 'bwr')
plt.colorbar()
plt.show()
