import os
from glob import glob

data_folder = "/mnt/grey/DATA/users/JB/test_HiM/test_data_Marion"
dest_folder = "/mnt/grey/DATA/users/JB/test_HiM/test_data_Marion/Deconvolved"


# to avoid pop up warning
huOpt.verb(mode = "silent")

# look for all the tif files
files = glob(data_folder + '/**/*.tif', recursive=True)
print(files)

# raw = image(path=)
print("DONE!")