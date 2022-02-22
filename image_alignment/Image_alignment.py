# -*- coding: utf-8 -*-
"""
@authors: JB. Fiche
Created on Mon January 10, 2022
------------------------------------------------------------------------------------------------------------------------
This code is used to align two DAPI images (image registration), taken at two different time points of an experiment,
after disassembly / reassembly of the fluidics chamber. The mathematical process used below is image correlation.
------------------------------------------------------------------------------------------------------------------------
"""

from tifffile import imread, TiffWriter
from scipy.signal import correlate
from scipy.ndimage import rotate, gaussian_filter, shift
# from time import time
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
# from skimage.filters import threshold_yen, threshold_multiotsu
# from skimage.exposure import rescale_intensity

matplotlib.use('TkAgg')


class ImageAlignment:

    def __init__(self, downsizing_power=3, gaussian_filter=10, angles_range=(-90, 90), crude_search_angle=9,
                 verbose=False, saving=False):
        self.downsizing_power = downsizing_power
        self.gaussian_filter = gaussian_filter
        self.angles_range = angles_range
        self.crude_search_angle = crude_search_angle
        self.verbose = verbose
        self.saving = saving
        self.resize_shape = None
        self.saving_path = ""
        self.optimum_angle = None
        self.dx = None
        self.dy = None

    def load_images(self, path_1, path_2, saving_path):
        """ Load the two images that are going to be realigned and calculate the parameters.

        :param saving_path: path to the folder where the results will be saved
        :param path_1: path leading to the first image (use as a reference)
        :param path_2: path leading to the second image that will be realigned with respect to image_1
        :return: im1, im2 the two loaded images.
        """
        im1 = imread(path_1)
        mip_1 = np.max(im1, axis=0)
        mip_1 = mip_1.astype(np.float_)

        im2 = imread(path_2)
        mip_2 = np.max(im2, axis=0)
        mip_2 = mip_2.astype(np.float_)

        self.saving_path = saving_path

        return mip_1, mip_2

    def process_image(self, im, im_name=None):
        """ Calculate the standardized image applying binning. The binning (down-sizing) is used to keep the essential
        structural information while saving time later during the alignment process. The correlation is indeed much
        faster when working on small images.

        :param im_name: name of the file in case the saving option is selected
        :param im: input image stack
        :return: mip : maximum intensity projection of the original stack im
                 im_standardized : return the standardized mip after binning
        """

        # Bin the image according to the downsizing power indicated as parameter
        self.resize_shape = (
                        int(im.shape[0] / 2 ** self.downsizing_power), int(im.shape[1] / 2 ** self.downsizing_power))

        shape = (self.resize_shape[0], im.shape[0] // self.resize_shape[0],
                 self.resize_shape[1], im.shape[1] // self.resize_shape[1])
        mip_bin = im.reshape(shape).mean(-1).mean(1)

        # Correct for illumination inhomogeneity using a gaussian filter
        if self.downsizing_power > 0:
            filter_size = self.gaussian_filter / (2 * self.downsizing_power)
            im_processed = np.divide(mip_bin, gaussian_filter(mip_bin, filter_size))
        else:
            im_processed = np.divide(mip_bin, gaussian_filter(mip_bin, self.gaussian_filter))

        # Standardized the image
        im_standardized = (im_processed - np.mean(im_processed)) / np.std(im_processed)

        # Save the images if the option was selected
        if self.saving and im_name:
            mip_name = im_name + '_MIP.tif'
            self.im_save(im, mip_name)

            processed_name = im_name + '_PROCESSED.tif'
            self.im_save(im_processed, processed_name)

        return im_standardized

    def alignment(self, im_ref, im, crude_search=True):
        """ Perform a search procedure to align an input image im with respect to im_ref. The aim is to find the
        rotation angle leading to the highest correlation score for the two images. The procedure can work two ways :
        either in a crude mode, where a crude search is performed according to the specified parameters. Or a fine
        search is run (with an angle step of 1°).

        :param im_ref: reference 2D input image
        :param im: image to be aligned
        :param crude_search: if True, a crude search is run
        :return: the optimum angle value (returning the highest correlation value) as well as the optimum shift values
        """
        # Select the central roi of the image that needs to be realigned
        roi, x0, y0 = self.select_roi(im)
        roi_shape = roi.shape[0]

        # Define the list of angle values to test
        if crude_search:
            n_angles = (self.angles_range[1] - self.angles_range[0]) // self.crude_search_angle
            angles = np.linspace(self.angles_range[0], self.angles_range[1], num=n_angles, endpoint=False,
                                 retstep=False,
                                 dtype=None, axis=0)
        elif crude_search is False and not self.optimum_angle:
            n_angles = self.angles_range[1] - self.angles_range[0] + 1
            angles = np.linspace(self.angles_range[0], self.angles_range[1], num=n_angles, endpoint=False,
                                 retstep=False,
                                 dtype=None, axis=0)
        else:
            angles = np.arange(self.optimum_angle - (self.crude_search_angle - 1),
                               self.optimum_angle + (self.crude_search_angle - 1), 1)
            n_angles = len(angles)

        correlation_max = np.zeros((n_angles,))
        correlation_im = np.zeros((n_angles, im_ref.shape[0], im_ref.shape[1]))
        max_correlation_2d_pos = np.zeros((n_angles, 2))

        # For each angle, rotate the roi and calculate the correlated image with the reference
        for n, angle in enumerate(angles):
            roi_rotated = self.im_rotate(roi, roi_shape, angle)
            corr = correlate(im_ref, roi_rotated, mode='same')
            correlation_max[n] = np.max(corr)
            correlation_im[n] = corr
            max_correlation_2d_pos[n, :] = np.unravel_index(np.argmax(corr, axis=None), corr.shape)

        # Plot the correlated images if the verbose option was selected
        if self.verbose:
            self.plot_correlation(correlation_im, angles)

        # Return the angle with the best correlation score as well as the shift dx, dy values
        max_correlation_idx = np.argmax(correlation_max)
        self.optimum_angle = angles[max_correlation_idx]
        self.dx = (x0 - max_correlation_2d_pos[max_correlation_idx, 0]) * 2 ** self.downsizing_power
        self.dy = (y0 - max_correlation_2d_pos[max_correlation_idx, 1]) * 2 ** self.downsizing_power

    def select_roi(self, im):
        """ Select the central roi of the input image.

        :param im: input 2D image
        :return: roi : the central part of the input image
        """
        d = int(self.resize_shape[0] / 2)
        x0 = int(self.resize_shape[0] / 2 - 1)
        y0 = int(self.resize_shape[1] / 2 - 1)
        roi = im[x0 - d // 2:x0 + d // 2, y0 - d // 2:y0 + d // 2]
        roi = (roi - np.mean(roi)) / np.std(roi)

        return roi, x0, y0

    @staticmethod
    def im_rotate(im, im_size, angle_value):
        """ Calculate the rotated image and crop the central part that is not affected by padding.

        :param im_size: size of the input image (assume to have the same x,y dimensions)
        :param im: input 2D image
        :param angle_value: value of the rotation
        :return: im_rotated, the cropped central part of the rotated image
        """

        # apply rotation by angle_value while keeping the same shape of the image
        im_rotated = rotate(im, angle_value, axes=(1, 0), reshape=False, order=3, mode='constant', prefilter=True)

        # select the central part of the rotated image (the one that will never be modified by padding whatever the
        # value of theta)
        d_shift = np.ceil((im_size - im_size / np.sqrt(2)) / 2)
        d_shift = d_shift.astype('int')
        im_rotated = im_rotated[d_shift:im_size - d_shift, d_shift:im_size - d_shift]

        return im_rotated

    def recalculate_image(self, im, processed=False):
        """ Apply the alignment transformation to the input 2D image.

        :param processed: indicate whether the input image has been processed (since it will change the dimensions of
        the image according to the downsizing parameter)
        :param im: 2D image (np array)
        :return: im_rotated_shifted, the modified image.
        """

        im = self.rescale_contrast(im)
        im_rotated = rotate(im, self.optimum_angle, axes=(1, 0), reshape=False, order=3, mode='constant', cval=0.0)
        if processed:
            im_rotated_shifted = shift(im_rotated,
                                       [-self.dx // 2 ** self.downsizing_power, -self.dy // 2 ** self.downsizing_power],
                                       order=0, mode='constant')
        else:
            im_rotated_shifted = shift(im_rotated, [-self.dx, -self.dy], order=0, mode='constant')

        return im_rotated_shifted

    def im_save(self, im, name):
        """ Save input image as a tiff file.

        :param im: input 2D images
        :param name: saving name
        """
        im_path = os.path.join(self.saving_path, name)
        im_min = np.min(im)
        im_max = np.max(im)
        im = (2 ** 16 - 1) * (im - im_min) / (im_max - im_min)
        with TiffWriter(im_path) as tif:
            tif.save(im.astype(np.uint16))

    @staticmethod
    def plot_correlation(im, angles):
        """ Plot all the correlated images if the verbose option is selected.

        :param im: array of images
        :param angles: array of angles tested for the correlation
        """
        plt.figure(figsize=(16, 10))
        num_im = im.shape[0]

        subplot_x = np.round(np.sqrt(num_im))
        subplot_y = np.ceil(num_im / subplot_x)

        for n in range(num_im):
            plt.subplot(subplot_x, subplot_y, n + 1)
            plt.imshow(im[n], cmap='gray')
            plt.title(f'{angles[n]}°')
            plt.axis('off')

        plt.show()

    def rescale_contrast(self, im):
        """ This function is used to improve the contrast of an input image, in order to make the visualization easier.

        :param im: input 2D image (np array)
        :return: 2D image with rescaled intensity
        """
        im = np.divide(im, gaussian_filter(im, self.gaussian_filter))
        intensity = np.copy(im)
        intensity = np.reshape(intensity, (intensity.shape[0] * intensity.shape[1], 1))
        intensity = np.sort(intensity, axis=0)
        int_min = np.percentile(intensity, 0.5)
        int_max = np.percentile(intensity, 99.5)
        im = (2 ** 16 - 1) * (im - int_min) / (int_max - int_min)
        im[im < 0] = 0
        im[im > 2**16-1] = 2**16
        return im

    def save_aligned_montage(self, im_ref, im_2_align, im_aligned, im_name="MIP"):
        """ Save the results.y

        :param im_ref: reference 2D image
        :param im_2_align: image to align
        :param im_aligned: aligned image
        :param im_name: name of the image to save
        """
        plt.figure(figsize=(16, 10))

        # Save the images as a subplot
        plt.subplot(131)
        plt.imshow(im_ref, cmap='gray')
        plt.title('Reference')
        plt.axis('off')
        plt.subplot(132)
        plt.imshow(im_2_align, cmap='gray')
        plt.title('Original image to align')
        plt.axis('off')
        plt.subplot(133)
        plt.imshow(im_ref, cmap='Blues')
        plt.imshow(im_aligned, cmap='Reds', alpha=0.5)
        plt.title(f'Aligned images - angle {self.optimum_angle}° - shift {self.dx, self.dy}px' )
        plt.axis('off')

        im_title = im_name + '_aligned.png'
        im_path = os.path.join(self.saving_path, im_title)
        plt.savefig(im_path)

        # Save the aligned montage :
        im_title = im_name + '_aligned_montage.png'
        im_path = os.path.join(self.saving_path, im_title)
        plt.figure(figsize=(16, 10))
        plt.imshow(im_ref, cmap='Blues')
        plt.imshow(im_aligned, cmap='Reds', alpha=0.5)
        plt.title(f'Aligned images - angle {self.optimum_angle}° - shift {self.dx, self.dy}px' )
        plt.axis('off')
        plt.savefig(im_path)


if __name__ == "__main__":
    # Indicate the path to the two images. The first image would be the reference. The second image is going to be
    # realigned with respect to the first.
    # path_image_ref = "/home/jb/Desktop/HiM_alignment/Test_2/HiM_DAPI/ROI_006/DAPI/scan_001_DAPI_006_ROI.tif"
    # path_image_to_align = "/home/jb/Desktop/HiM_alignment/Test_2/RNA_DAPI/DAPI_raw/ROI_006/DAPI/scan_001_DAPI_006_ROI.tif"
    # dest_folder = "/home/jb/Desktop/HiM_alignment/Test_2/ROI_6_results"

    # path_image_ref = "/home/jb/Desktop/HiM_alignment/Test_2/HiM_DAPI/ROI_015/DAPI/scan_001_DAPI_015_ROI.tif"
    # path_image_to_align = "/home/jb/Desktop/HiM_alignment/Test_2/RNA_DAPI/DAPI_raw/ROI_015/DAPI/scan_001_DAPI_015_ROI.tif"
    # dest_folder = "/home/jb/Desktop/HiM_alignment/Test_2/ROI_15_results"

    path_image_ref = "/home/jb/Desktop/HiM_alignment/Test_1/DAPI_OM_MS_1.tif"
    path_image_to_align = "/home/jb/Desktop/HiM_alignment/Test_1/DAPI_OM_MS_2.tif"
    dest_folder = "/home/jb/Desktop/HiM_alignment/Test_1"

    # Instantiate the alignment class
    _align = ImageAlignment(downsizing_power=3, gaussian_filter=15, angles_range=(-90, 90), crude_search_angle=9,
                            verbose=False, saving=True)

    # Load the images and define the main parameters
    mip_reference, mip_to_align = _align.load_images(path_image_ref, path_image_to_align, dest_folder)

    # Process the two images
    im_reference_processed = _align.process_image(mip_reference, im_name="im_ref")
    im_to_align_processed = _align.process_image(mip_to_align, im_name="im_to_align")

    # # Perform a first crude search for the alignment
    # _align.alignment(im_reference_processed, im_to_align_processed, crude_search=True)
    #
    # # Based on the previous calculation, perform a second alignment with constrained angle values
    # _align.alignment(im_reference_processed, im_to_align_processed, crude_search=False)
    #
    # # Calculate the final image using the MIP
    # aligned_image = _align.recalculate_image(mip_to_align, processed=False)
    # mip_ref_contrast = _align.rescale_contrast(mip_reference)
    # mip_to_align_contrast = _align.rescale_contrast(mip_to_align)
    # _align.save_aligned_montage(mip_ref_contrast, mip_to_align_contrast, aligned_image, im_name="MIP")
