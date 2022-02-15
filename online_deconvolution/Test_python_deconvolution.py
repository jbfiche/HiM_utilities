from glob import glob
import os
import sys
version = f'Python version : {sys.version}'
huOpt.report(version)

# parameters
data_folder = "/mnt/grey/DATA/users/JB/test_HiM/test_data_Marion/Raw_data"
dest_folder = "/mnt/grey/DATA/users/JB/test_HiM/test_data_Marion/Deconvolved"

n_channel: int = 2
wavelength: list = [561, 647]
NA: float = 1.2
voxel_size: list = [106, 106, 250] #nm

# to avoid pop up warning
huOpt.verb(mode="silent")

# define the psf according to the properties define in parameters


# look for all the tif files
path_files = glob(data_folder + '/**/*.tif', recursive=True)
text_report = f'Number of tif files found : {len(path_files)}'
huOpt.report(text_report)

# for each tif file, open an image class
for n, path in enumerate(path_files):
    text_report = f'Loading image #{n} with the following path : {path}'
    huOpt.report(text_report)
    raw = image(path=path, foreignTo="float", series="off")
    
    # get dimensions of the image
    im_size = raw.getdims(mode="all")
    dX = im_size[0]
    dY = im_size[1]
    n_frames = im_size[2]
    n_ch_frames = int(n_frames / n_channel)
    text_report = f'Image dimensions : {im_size} and number of frames : {n_frames}'
    huOpt.report(text_report)
    
    # depending on the number of channels, split the data accordingly
    for ch in range(n_channel):
        raw_channel = image("img_channel", type = "float", dim = [dX, dY, n_ch_frames, 0, 1, 1])
        
        # copy from raw the frames associated to the selected channel
        for frame in range(n_ch_frames):
            src_frame = int(n_channel*frame + ch)
            raw.cp(raw_channel, span = [dX, dY, 1, 0, 1, 1], srco = [0, 0, src_frame, 0, 0, 0], desto = [0, 0, frame, 0, 0, 0])

        # perform deconvolution
        
        # save the deconvolved image
        file_name = os.path.basename(path)
        new_file_name = os.path.splitext(file_name)[0] + "_decon_ch" + str(ch) + ".tif"
        saving_path = os.path.join(dest_folder, new_file_name)
        raw_channel.save(saving_path, type="tiff16", tiffMultiDir=True, cmode="scale")

huOpt.report("Done!")
exit()