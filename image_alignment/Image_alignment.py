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
from time import time
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt

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
        self.im_shape = None
        self.resize_shape = None
        self.saving_path = ""

    def load_images(self, path_1, path_2, saving_path):
        """ Load the two images that are going to be realigned and calculate the parameters.

        :param saving_path: path to the folder where the results will be saved
        :param path_1: path leading to the first image (use as a reference)
        :param path_2: path leading to the second image that will be realigned with respect to image_1
        :return: im1, im2 the two loaded images.
        """
        im1 = imread(path_1)
        im2 = imread(path_2)
        self.im_shape = im1.shape
        self.saving_path = saving_path

        return im1, im2

    def process_image(self, im, im_name=None):
        """ Calculate the standardized image applying binning. The binning is applied after calculating the MIP in order to
    optimize the calculation time. The image renormalization using the gaussian filter is actually quite long when working
    on full size images

        :param im_name: name of the file in case the saving option is selected
        :param im: input image stack
        :return: mip : maximum intensity projection of the original stack im
                 im_standardized : return the standardized mip after binning
        """
        # Calculate the mip
        mip = np.max(im, axis=0)
        mip = mip.astype(np.float_)

        # Bin the image according to the downsizing power indicated as parameter
        self.resize_shape = (int(self.im_shape[1] / 2 ** self.downsizing_power), int(self.im_shape[2] / 2 ** self.downsizing_power))
        shape = (self.resize_shape[0], mip.shape[0] // self.resize_shape[0],
                 self.resize_shape[1], mip.shape[1] // self.resize_shape[1])
        mip_bin = mip.reshape(shape).mean(-1).mean(1)

        # Correct for illumination inhomogeneity using a gaussian filter
        if self.downsizing_power > 0:
            self.gaussian_filter = 10 / (2 * self.downsizing_power)
        im_processed = np.divide(mip_bin, gaussian_filter(mip_bin, self.gaussian_filter))

        # Standardized the image
        im_standardized = (im_processed - np.mean(im_processed)) / np.std(im_processed)

        # Save the images if the option was selected
        if self.saving and im_name:
            mip_name = im_name + '_MIP.tif'
            self.im_save(mip, mip_name)

            processed_name = im_name + '_PROCESSED.tif'
            self.im_save(im_processed, processed_name)

        return im_standardized

    def alignment(self, im_ref, im, crude_search=True, starting_angle=None):
        """ Perform a search procedure to align an input image im with respect to im_ref. The aim is to find the
        rotation angle leading to the highest correlation score for the two images. The procedure can work two ways :
        either in a crude mode, where a crude search is performed according to the specified parameters. Or a fine
        search is run (with an angle step of 1°).

        :param im_ref: reference 2D input image
        :param im: image to be aligned
        :param crude_search: if True, a crude search is run
        :param starting_angle: if indicated, the fine search will be performed around this values
        :return: the optimum angle value (returning the highest correlation value)
        """
        # Select the central roi of the image that needs to be realigned
        roi = self.select_roi(im)
        roi_shape = roi.shape[0]

        # Define the list of angle values to test
        if crude_search:
            n_angles = (self.angles_range[1] - self.angles_range[0]) // self.crude_search_angle
            angles = np.linspace(self.angles_range[0], self.angles_range[1], num=n_angles, endpoint=False, retstep=False,
                                 dtype=None, axis=0)
        elif crude_search is False and not starting_angle:
            n_angles = self.angles_range[1] - self.angles_range[0] + 1
            angles = np.linspace(self.angles_range[0], self.angles_range[1], num=n_angles, endpoint=False, retstep=False,
                                 dtype=None, axis=0)
        else:
            angles = np.arange(starting_angle - (self.crude_search_angle - 1),
                               starting_angle + (self.crude_search_angle - 1), 1)
            n_angles = len(angles)

        correlation_max = np.zeros((n_angles,))
        correlation_im = np.zeros((n_angles, im_ref.shape[0], im_ref.shape[1]))

        # For each angle, rotate the roi and calculate the correlated image with the reference
        for n, angle in enumerate(angles):
            roi_rotated = self.im_rotate(roi, roi_shape, angle)
            corr = correlate(im_ref, roi_rotated, mode='same')
            correlation_max[n] = np.max(corr)
            correlation_im[n] = corr

        # Plot the correlated images if the verbose option was selected
        if self.verbose:
            self.plot_correlation(correlation_im, angles)

        # Return the angle with the best correlation score
        return angles[np.argmax(correlation_max)]

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

        return roi

    @staticmethod
    def im_rotate(im, im_size, angle_value):
        """ Calculate the rotated image and crop the central part that is not affected by padding.

        :param im_size: size of the input image (assume to have the same x,y dimensions)
        :param im: input 2D image
        :param angle_value: value of the rotation
        :return: im_rotated, the cropped central part of the rotated image
        """

        # apply rotation by Theta while keeping the same shape of the image
        im_rotated = rotate(im, angle_value, axes=(1, 0), reshape=False, order=3, mode='constant', prefilter=True)

        # select the central part of the rotated image (the one that will never be modified by padding whatever the
        # value of theta)
        d_shift = np.ceil((im_size - im_size / np.sqrt(2)) / 2)
        d_shift = d_shift.astype('int')
        im_rotated = im_rotated[d_shift:im_size - d_shift, d_shift:im_size - d_shift]

        return im_rotated

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


# # Select the central ROI from DAPI_1 :
# d = roi_size
# x0 = int(resize_shape[0] / 2 - 1)
# y0 = int(resize_shape[0] / 2 - 1)
# roi = DAPI_1_processed[x0 - d // 2:x0 + d // 2, y0 - d // 2:y0 + d // 2]
# roi = (roi - np.mean(roi)) / np.std(roi)
# DAPI_1_standardize = (DAPI_1_processed - np.mean(DAPI_1_processed)) / np.std(DAPI_1_processed)
# DAPI_2_standardize = (DAPI_2_processed - np.mean(DAPI_2_processed)) / np.std(DAPI_2_processed)


# # Calculate the correlation of the roi with DAPI_2, after applying a rotation defined by the angle theta. Crude search
# # with a precision of 10° :
# def im_correlate(im1, im2, angle_value, im_size):
#     # apply rotation by Theta while keeping the same shape of the image
#     im1_rotate = rotate(im1, angle_value, axes=(1, 0), reshape=False, order=3, mode='constant', prefilter=True)
#     # select the central part of the image, the one that will never be modified by padding whatever the value of theta
#     d_shift = np.ceil((im_size - im_size / np.sqrt(2)) / 2)
#     d_shift = d_shift.astype('int')
#     im1_rotate = im1_rotate[d_shift:im_size - d_shift, d_shift:im_size - d_shift]
#     # calculate the correlation between the two images
#     corr = correlate(im2, im1_rotate, mode='same')
#     return corr



# n_angles = (angles_range[1] - angles_range[0]) // crude_search_angle
# angles = np.linspace(angles_range[0], angles_range[1], num=n_angles, endpoint=False, retstep=False, dtype=None, axis=0)
# correlation_max = np.zeros((360 // crude_search_angle,))
# if verbose:
#     plt.figure(figsize=(16, 10))
#
# for n, angle in enumerate(angles):
#     correlation_image = im_correlate(roi, DAPI_2_standardize, angle, d)
#     correlation_max[n] = np.max(correlation_image)
#     # In case the option is selected, display all the images in one single figure
#     if verbose:
#         im_plot(correlation_image, n, angle)
#
# if verbose:
#     plt.show()

# # From the maximum correlation found in the first crude search, perform a finer search with a precision of 1°. Infer the
# # best angle for the image rotation and calculate the shift for the translation :
# theta_start = angles[np.argmax(correlation_max)] - (crude_search_angle - 1)
# theta_stop = angles[np.argmax(correlation_max)] + (crude_search_angle - 1)
# angles = np.arange(theta_start, theta_stop, 1)
# correlation_max = np.zeros((len(angles),))
# index_max = np.zeros((len(angles), 2))
#
# for n, angle in enumerate(angles):
#     correlation_image = im_correlate(roi, DAPI_2_standardize, angle, d)
#     correlation_max[n] = np.max(correlation_image)
#     index_max[n, :] = np.unravel_index(np.argmax(correlation_image, axis=None), correlation_image.shape)
#
# theta = angles[np.argmax(correlation_max)]
# dx = x0 - index_max[np.argmax(correlation_max), 0]
# dy = y0 - index_max[np.argmax(correlation_max), 1]
#
# t3 = time()
# print(f'Optimization time : {t3 - t2}s')
#
# # Calculate the final image and the total calculation time :
# DAPI_1_rotate = rotate(DAPI_1_mip, theta, axes=(1, 0), reshape=False, order=3, mode='constant', cval=0.0)
# DAPI_1_rotate_shift = shift(DAPI_1_rotate, [-dx * 2 ** downsize_power, -dy * 2 ** downsize_power], order=0,
#                             mode='constant')
#
# t4 = time()
# print(f'Registration of the original image : {t4 - t3}s')
# print(f'Total calculation time : {t4 - t0}s')
#
# # Plot the final images using the MIP:
# plt.figure(figsize=(16, 10))
# plt.subplot(131)
# plt.imshow(DAPI_1_mip, cmap='gray')
# plt.title('DAPI_1')
# plt.axis('off')
# plt.subplot(132)
# plt.imshow(DAPI_2_mip, cmap='gray')
# plt.title('DAPI_2')
# plt.axis('off')
# plt.subplot(133)
# plt.imshow(DAPI_2_mip, cmap='gray')
# plt.imshow(DAPI_1_rotate_shift, cmap='twilight', alpha=0.5)
# plt.title('Aligned images (DAPI_2 as reference)')
# plt.axis('off')
#
# im_name = os.path.join(current_path, 'MIP_aligned.png')
# plt.savefig(im_name)
#
# # Plot the final images using the standardized images:
# DAPI_1_rotate = rotate(DAPI_1_processed, theta, axes=(1, 0), reshape=False, order=3, mode='constant', cval=0.0)
# DAPI_1_rotate_shift = shift(DAPI_1_rotate, [-dx, -dy], order=3, mode='constant')
#
# plt.figure(figsize=(16, 10))
# plt.subplot(131)
# plt.imshow(DAPI_1_processed, cmap='gray')
# plt.title('DAPI_1')
# plt.axis('off')
# plt.subplot(132)
# plt.imshow(DAPI_2_processed, cmap='gray')
# plt.title('DAPI_2')
# plt.axis('off')
# plt.subplot(133)
# plt.imshow(DAPI_2_processed, cmap='gray')
# plt.imshow(DAPI_1_rotate_shift, cmap='twilight', alpha=0.5)
# plt.title('Aligned images (DAPI_2 as reference - DAPI_1 in blue)')
# plt.axis('off')
#
# im_name = os.path.join(current_path, 'MIP_corrected_aligned.png')
# plt.savefig(im_name)
#
# # Save the aligned montage :
# im_name = os.path.join(current_path, 'Aligned_images.png')
# plt.figure(figsize=(16, 10))
# plt.imshow(DAPI_2_processed, cmap='gray')
# plt.imshow(DAPI_1_rotate_shift, cmap='twilight', alpha=0.5)
# plt.title('Aligned images (DAPI_2 as reference - DAPI_1 in blue)')
# plt.axis('off')
# plt.savefig(im_name)

if __name__ == "__main__":

    # Indicate the path to the two images. The first image would be the reference. The second image is going to be
    # realigned with respect to the first.
    path_image_1 = "/home/jb/Desktop/HiM_alignment/Test_2/HiM_DAPI/ROI_006/DAPI/scan_001_DAPI_006_ROI.tif"
    path_image_2 = "/home/jb/Desktop/HiM_alignment/Test_2/RNA_DAPI/DAPI_raw/ROI_006/DAPI/scan_001_DAPI_006_ROI.tif"
    dest_folder = "/home/jb/Desktop/HiM_alignment/Test_2"

    # Instantiate the alignment class
    _align = ImageAlignment(downsizing_power=3, gaussian_filter=10, angles_range=(-90, 90), crude_search_angle=9,
                            verbose=True, saving=True)

    # Load the images and define the main parameters
    im_reference, im_to_align = _align.load_images(path_image_1, path_image_2, dest_folder)

    # Process the two images
    im_reference_processed = _align.process_image(im_reference, im_name="im_ref")
    im_to_align_processed = _align.process_image(im_to_align, im_name="im_to_align")

    # Perform a first crude search for the alignment
    best_angle_approximation = _align.alignment(im_reference_processed, im_to_align_processed, crude_search=True)

    # Based on the previous calculation, perform a second alignment with constrained angle values
    best_angle = _align.alignment(im_reference_processed, im_to_align_processed, crude_search=False,
                                  starting_angle=best_angle_approximation)



