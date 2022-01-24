import os
from glob import glob

data_folder = "/mnt/grey/DATA/users/JB/test_HiM/test_data_Marion"
dest_folder = "/mnt/grey/DATA/users/JB/test_HiM/test_data_Marion/Deconvolved"


# to avoid pop up warning
huOpt.verb(mode="silent")

# look for all the tif files
path_files = glob(data_folder + '/**/*.tif', recursive=True)
print(f'Number of tif files found : {len(path_files)}')

# for each tif file, open an image class
for n, path in enumerate(path_files):
    print(f'Loading image #{n} with the following path : {path}')
    raw = image(path=path)

print("DONE!")