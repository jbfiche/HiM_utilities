#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 17:43:50 2021

@author: fiche
"""

import numpy as np
from glob import glob
from tifffile import imread, TiffWriter
import os
import time
#import matplotlib.pyplot as plt

# data_folder = "/mnt/grey/DATA/rawData_2020/Experiment_18_Olivier_Locus2L_low_res/DAPI_raw/002_merFISH_RAMM_DAPI"
#data_folder = "/mnt/grey/DATA/rawData_2021/Experiment_40_Marion/Brain/2021_08_05/003_HiM_Hires_lib_brain"
data_folder = "/mnt/grey/DATA/2021_09_17/001_HiM_Tissues_Olivier_HiM/"
dest_folder = "/mnt/grey/DATA/2021_09_17/Deconvolved/"

list_analyzed = []
last_deconvolution_time = time.time()

# Launch a while loop. This loop stops if no new data is found for 3h
# -------------------------------------------------------------------

while time.time() - last_deconvolution_time < 10800:
    
    list_to_analyze = []

# List all the .tif files in the selected folder and sort them according to 
# their creation time. By security, we check that the size of the most recent
# file is not changing after a couple of second. 
# ----------------------------------------------

    os.chdir(data_folder)
    im_data = glob('**/*.tif', recursive=True)
    im_data.sort(key=os.path.getctime)
    
    if im_data:
        size_t0 = os.path.getsize(im_data[-1])
        time.sleep(5)
        size_t1 = os.path.getsize(im_data[-1])
        if size_t0 != size_t1:
            del im_data[-1]
        
    print(im_data)
    
    # Compare the list im_data with the list containing the names of the previously 
    # analyzed images
    # ---------------
    
    for im in im_data:
        if im not in list_analyzed:
            list_to_analyze.append(im)
    
    if list_to_analyze :
        
        print("List of the images to be processed : {}".format(list_to_analyze))
    
        # Modify the tcl script for the deconvolution in order to perform the operation 
        # on a list of images
        # -------------------
        
        os.chdir('/home/fiche/huygensDeconvolutionScripts/')
        singlefile_decon_time = np.zeros(len(im_data),)
        
        a_file = open("batchDeconvolve_online.tcl", "r")
        list_of_lines = a_file.readlines()
        list_of_lines[5] = "set destDir " + dest_folder + " \n"
        
        im_list = "set dataList [list "
        
        for n,im in enumerate(list_to_analyze):
            im_list = im_list + os.path.join(data_folder, im) + " "
            
        im_list = im_list + "] \n"
        list_of_lines[3] = im_list
        
        a_file = open("batchDeconvolve_online.tcl", "w")
        a_file.writelines(list_of_lines)
        a_file.close()
        
        # Launch the deconvolution
        # ------------------------
        
        start = time.time()
        # os.system('nohup hucore -task /home/fiche/huygensDeconvolutionScripts/batchDeconvolve_online.tcl > out.log')
        os.system('hucore -task /home/fiche/huygensDeconvolutionScripts/batchDeconvolve_online.tcl')
        end = time.time()
        print((end - start)/60)
        
        last_deconvolution_time = end
        list_analyzed = list_to_analyze + list_analyzed
        print("List of the processed images : {}".format(list_analyzed))
        
    else :
        
        dt = ( time.time() - last_deconvolution_time ) / 60
        print('No new file was found since {} min'.format(np.round(dt)))
        time.sleep(5)



# for n,im in enumerate(im_data):
#     a_file = open("batchDeconvolve_online.tcl", "r")
#     list_of_lines = a_file.readlines()
    
#     list_of_lines[3] = "set dataDir " + os.path.join(data_folder, os.path.split(im)[0]) + "/ \n"
#     list_of_lines[5] = "set destDir " + dest_folder + " \n"
     
#     a_file = open("batchDeconvolve.tcl", "w")
#     a_file.writelines(list_of_lines)
#     a_file.close()
    
#     start = time.time()
#     # os.system('nohup hucore -task /home/fiche/huygensDeconvolutionScripts/batchDeconvolve.tcl > out.log')
#     os.system('hucore -task /home/fiche/huygensDeconvolutionScripts/batchDeconvolve.tcl')
#     end = time.time()
#     print((end - start)/60)
#     singlefile_decon_time[n] = (end - start)/60
    
# # Plot the analysis time
# # ----------------------

# fig = plt.figure(figsize=(8,5))
# plt.plot(singlefile_decon_time,   'o-')
# plt.xlabel('Iterations')
# plt.ylabel('Deconvolution time (min)')

# fig_name = dest_folder + '/' + 'Deconvolution_time.png'
# fig.savefig(fig_name)
