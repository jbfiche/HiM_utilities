from glob import glob
import os
import sys
version = f'Python version : {sys.version}'
huOpt.report(version)

# parameters
data_folder = "/mnt/grey/DATA/users/JB/test_HiM/Test_data_Olivier_embryo/raw_DAPI"
dest_folder = "/mnt/grey/DATA/users/JB/test_HiM/Test_data_Olivier_embryo/Results_python"

n_channel: int = 2
wavelength_excitation: list = [405, 561, 670]
wavelength_emission: list = [460, 575, 700]  # Emission maxima for DAPI, Atto550, AF647
NA: float = 1.2
voxel_size: list = [106, 106, 250]  # nm
ri_sample: float = 1.33  # Specimen refractive index
ri_medium: float = 1.515  # Lens medium refractive index
objQuality: str = "good"  # Quality of the objective lens
microscope: str = "widefield"  # Imaging modality
verbose: bool = True

# to avoid pop up warning
huOpt.verb(mode="silent")

# define the psf according to the properties define in parameters
psf = image("psf", type="float")
psf = psf.genpsfExpl(na=NA, ri=ri_sample, ril=ri_medium, ex=wavelength_excitation[0],
                     em=wavelength_emission[0], dims="auto", dx=50, dz=100, micr=microscope, zDist=0.0,
                     imagingDir="upward", reflCorr=False, objQuality=objQuality, v=verbose)

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
        raw_channel = image("img_channel", type="float", dim=[dX, dY, n_ch_frames, 0, 1, 1])

        # copy from raw the frames associated to the selected channel
        for frame in range(n_ch_frames):
            src_frame = int(n_channel*frame + ch)
            raw.cp(raw_channel, span=[dX, dY, 1, 0, 1, 1], srco=[0, 0, src_frame, 0, 0, 0], desto=[0, 0, frame, 0, 0, 0])

        # perform deconvolution
        # im_param = raw_channel.setp(na=NA, objQuality=objQuality, ril=ri_medium, ri=ri_sample, ex=wavelength_excitation[0],
        #                  em=wavelength_emission[0], baseline=100, dx=0.106, dy=0.106, dz=0.25, mag=60.0,
        #                  imagingDir="upward", micr=microscope, tclReturn=True)
        # text_report = f'image paramters : {im_param}'

        deconvolved = raw_channel.cmle(psf, sn=[20, 20, 20, 20, 20], snr=[12, 12, 12, 12, 12], it=40, bgMode="wf",
                                       bg=[0.0, 0.0, 0.0, 0.0], blMode="off", brMode="auto", varPsf="off", q=0.1,
                                       mode="fast", pad="auto", reduceMode="auto", bgRadius=0.7)

        # save the deconvolved image
        file_name = os.path.basename(path)
        new_file_name = os.path.splitext(file_name)[0] + "_decon_ch" + str(ch) + ".tif"
        saving_path = os.path.join(dest_folder, new_file_name)
        deconvolved.save(saving_path, type="tiff16", tiffMultiDir=True, cmode="scale")

huOpt.report("Done!")
exit()