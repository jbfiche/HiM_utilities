#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 16:23:45 2020

The following program is used to analyze the data from the Airyscan microscope,
rename the movies and convert them to tif files. In order to work, use the 
conda environment "Confocal"

@author: jb
"""

import os
import numpy as np
import re
import czifile
import tifffile

# import matplotlib.pyplot as plt

# Define the folders where the acquisition data were saved and where the 
# template txt file can be found
# -------------------------------

path_acquisition_data = ["/mnt/grey/DATA/rawData_2021/Experiment_64_Jb_Christophe_Pia/Frozen cells_2021-11-26/RT1/RT_control-1_AcquisitionBlock2.czi",
                         "/mnt/grey/DATA/rawData_2021/Experiment_64_Jb_Christophe_Pia/Frozen cells_2021-11-26/RT1/RT_control-1_AcquisitionBlock3.czi"]
path_name_file = "/mnt/grey/DATA/rawData_2021/Experiment_64_Jb_Christophe_Pia/Frozen cells_2021-11-26/rt1_movie_name.txt"
saving_path = "/mnt/grey/DATA/rawData_2021/Experiment_64_Jb_Christophe_Pia/Frozen cells_2021-11-26/RT1_"


# List the content of all the acquisition folders and save it in the list
# "Data_file". Since the name of the files depends on the prefix indicated by 
# the user, the acquisition point is retrieved by looking for the pattern 
# "_ptxx". Using the acqusition point, the files can be sorted according to 
# the acquisition order and the paths are saved in the list "Sorted_path". 
# ------------------------------------------------------------------------

Data_file = []
Data_number = []

for path in path_acquisition_data:
    for file in os.listdir(path):
        Data_file.append(os.path.join(path, file))
        n = re.findall(r'_pt\d+', file)
        Data_number.append(n[0][3:])

Data_number = np.array(Data_number)  
Data_number = Data_number.astype(np.int16)

Idx = np.argsort(Data_number)

Sorted_path = []
for n in Idx:
    Sorted_path.append(Data_file[n])
    
# Open the txt file and save the files names in the list "Data_name". Make
# sure the number of files in the acquisition folder match the number of items
# indicated in the txt file.
# --------------------------

with open(path_name_file) as f:
    Data_name = f.readlines()
    Data_name = [x.strip() for x in Data_name]
    
if len(Data_name) != len(Sorted_path):
    print(
        "There is an error in the selected folders. The number of acquisition movies does not match the number of file name in the txt file")

# Open the czi movies and save it as tiff

for nmovie in range(len(Sorted_path)):
    
    movie_name = Data_name[nmovie]
    delimiter_pos = []
    nchar = 0
    
    for char in movie_name:
        if char == "_":
           delimiter_pos.append(nchar)
        nchar = nchar+1
        
    prefix = movie_name[0:delimiter_pos[0]]
    embryo = movie_name[delimiter_pos[0]+1:delimiter_pos[1]]
    roi = movie_name[delimiter_pos[1]+1:]

    print(f'movie_name : {movie_name}, prefix : {prefix}, embryo : {embryo}, roi : {roi}')

    movie = czifile.imread(Sorted_path[nmovie])
    n_channel = movie.shape[1]
    n_z = movie.shape[2]
    
    tiff_name = 'Scan_' + embryo + '_' + roi + '.tif'
    tiff_path = os.path.join(saving_path, tiff_name)
    with tifffile.TiffWriter(tiff_path) as tf:
    
        for z in range(n_z):
            for c in range (n_channel):
        
                im = movie[0, c, z, :, :, 0]
                im = np.array(im)
                tf.save(im.astype(np.uint16))
            
    message = movie_name + " was saved"
    print(message)
