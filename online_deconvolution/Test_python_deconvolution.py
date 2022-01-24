import os
from glob import glob

data_folder = "/mnt/grey/DATA/users/JB/test_HiM/test_data_Marion/Test"
dest_folder = "/mnt/grey/DATA/users/JB/test_HiM/test_data_Marion/Deconvolved"

n_channel = 2

# to avoid pop up warning
huOpt.verb(mode="silent")

# look for all the tif files
path_files = glob(data_folder + '/**/*.tif', recursive=True)
# print(f'Number of tif files found : {len(path_files)} \n')
huOpt.report("Number of tif files found : ${len(path_files)}")

# create an empty image, so the deconvolution is done using a theoretical psf based on Microscopic Parameters
psf = image()

# for each tif file, open an image class
for n, path in enumerate(path_files):
    # print(f'Loading image #{n} with the following path : {path} \n')
    raw = image(path=path)
    result_img = raw.cmle(psf, it=5, bgMode="wf", bgRadius=0.7, blMode="off")

huOpt.report("Done!")
exit()