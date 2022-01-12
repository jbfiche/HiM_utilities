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

t0 = time()

# Define the main folder :
current_path = "/home/jb/Desktop"

# Load the data :
DAPI_1 = imread(os.path.join(current_path, "DAPI_OM_MS_1.tif"))
DAPI_2 = imread(os.path.join(current_path, "DAPI_OM_MS_2.tif"))
im_shape = DAPI_1.shape

t1 = time()
print(f'Loading images time : {t1-t0}s')

# Define the parameters
downsize_factor = 8  # Need to be a power of 2
gaussian_filter_size = 10 / np.log2(downsize_factor)
resize_shape = (int(im_shape[1] / downsize_factor), int(im_shape[2] / downsize_factor))
roi_size = int(resize_shape[0]/2)
angles_range = (-90, 90)
crude_search_angle = 9
verbose = False
if verbose:
    subplot_X = np.round(np.sqrt(360//crude_search_angle))
    subplot_Y = np.ceil((360//crude_search_angle)/subplot_X)
saving = False


# Calculate the standardized image applying binning. The binning is applied after calculating the MIP in order to
# optimize the calculation time. The image renormalization using the gaussian filter is actually quite long when working
# on full size images :
def standardize(im, new_shape):
    mip = np.max(im, axis=0)
    mip = mip.astype(np.float_)

    shape = (new_shape[0], mip.shape[0] // new_shape[0],
             new_shape[1], mip.shape[1] // new_shape[1])
    mip_bin = mip.reshape(shape).mean(-1).mean(1)

    im_processed = np.divide(mip_bin, gaussian_filter(mip_bin, gaussian_filter_size))
    im_min = np.min(im_processed)
    im_max = np.max(im_processed)
    im_standardized = (2**16 - 1) * (im_processed - im_min) / (im_max - im_min)
    return mip, im_standardized


DAPI_1_mip, DAPI_1_processed = standardize(DAPI_1, resize_shape)
DAPI_2_mip, DAPI_2_processed = standardize(DAPI_2, resize_shape)

t2 = time()
print(f'Standardization time : {t2-t1}s')


# Save the processed images :
def im_save(im, name, path):
    im_path = os.path.join(path, name)
    with TiffWriter(im_path) as tif:
        tif.save(im.astype(np.uint16))


if saving:
    im_save(DAPI_1_processed.astype(np.uint16), 'DAPI_1_processed', current_path)
    im_save(DAPI_2_processed.astype(np.uint16), 'DAPI_2_processed', current_path)

# Select the central ROI from DAPI_1 :
d = roi_size
x0 = int(resize_shape[0]/2 - 1)
y0 = int(resize_shape[0]/2 - 1)
roi = DAPI_1_processed[x0-d//2:x0+d//2, y0-d//2:y0+d//2]
roi = (roi - np.mean(roi))/np.std(roi)
DAPI_1_standardize = (DAPI_1_processed - np.mean(DAPI_1_processed))/np.std(DAPI_1_processed)
DAPI_2_standardize = (DAPI_2_processed - np.mean(DAPI_2_processed))/np.std(DAPI_2_processed)


# Calculate the correlation of the roi with DAPI_2, after applying a rotation defined by the angle theta. Crude search
# with a precision of 10° :
def im_correlate(im1, im2, angle_value, im_size):
    # apply rotation by Theta while keeping the same shape of the image
    im1_rotate = rotate(im1, angle_value, axes=(1, 0), reshape=False, order=3, mode='constant', prefilter=True)
    # select the central part of the image, the one that will never be modified by padding whatever the value of theta
    d_shift = np.ceil((im_size - im_size/np.sqrt(2))/2)
    d_shift = d_shift.astype('int')
    im1_rotate = im1_rotate[d_shift:im_size-d_shift, d_shift:im_size-d_shift]
    # calculate the correlation between the two images
    corr = correlate(im2, im1_rotate, mode='same')
    return corr


def im_plot(im, n_subplot, angle_value):
    plt.subplot(subplot_X, subplot_Y, n_subplot+1)
    plt.imshow(im, cmap='gray')
    plt.title(f'{angle_value}°')
    plt.axis('off')


n_angles = (angles_range[1] - angles_range[0]) // crude_search_angle
angles = np.linspace(angles_range[0], angles_range[1], num=n_angles, endpoint=False, retstep=False, dtype=None, axis=0)
correlation_max = np.zeros((360//crude_search_angle,))
if verbose:
    plt.figure(figsize=(16, 10))

for n, angle in enumerate(angles):
    correlation_image = im_correlate(roi, DAPI_2_standardize, angle, d)
    correlation_max[n] = np.max(correlation_image)
    # In case the option is selected, display all the images in one single figure
    if verbose:
        im_plot(correlation_image, n, angle)

if verbose:
    plt.show()

# From the maximum correlation found in the first crude search, perform a finer search with a precision of 1°. Infer the
# best angle for the image rotation and calculate the shift for the translation :
theta_start = angles[np.argmax(correlation_max)] - (crude_search_angle - 1)
theta_stop = angles[np.argmax(correlation_max)] + (crude_search_angle - 1)
angles = np.arange(theta_start, theta_stop, 1)
correlation_max = np.zeros((len(angles),))
index_max = np.zeros((len(angles), 2))

for n, angle in enumerate(angles):
    correlation_image = im_correlate(roi, DAPI_2_standardize, angle, d)
    correlation_max[n] = np.max(correlation_image)
    index_max[n, :] = np.unravel_index(np.argmax(correlation_image, axis=None), correlation_image.shape)

theta = angles[np.argmax(correlation_max)]
dx = x0 - index_max[np.argmax(correlation_max), 0]
dy = y0 - index_max[np.argmax(correlation_max), 1]

t3 = time()
print(f'Optimization time : {t3-t2}s')

# Calculate the final image and the total calculation time :
DAPI_1_rotate = rotate(DAPI_1_mip, theta, axes=(1, 0), reshape=False, order=3, mode='constant', cval=0.0)
DAPI_1_rotate_shift = shift(DAPI_1_rotate, [-dx * downsize_factor, -dy * downsize_factor], order=0, mode='constant')

t4 = time()
print(f'Registration of the original image : {t4-t3}s')
print(f'Total calculation time : {t4-t0}s')

# Plot the final images using the MIP:
plt.figure(figsize=(16, 10))
plt.subplot(131)
plt.imshow(DAPI_1_mip, cmap='gray')
plt.title('DAPI_1')
plt.axis('off')
plt.subplot(132)
plt.imshow(DAPI_2_mip, cmap='gray')
plt.title('DAPI_2')
plt.axis('off')
plt.subplot(133)
plt.imshow(DAPI_2_mip, cmap='gray')
plt.imshow(DAPI_1_rotate_shift, cmap='twilight', alpha=0.5)
plt.title('Aligned images (DAPI_2 as reference)')
plt.axis('off')
plt.show()

# Plot the final images using the standardized images:
DAPI_1_rotate = rotate(DAPI_1_processed, theta, axes=(1, 0), reshape=False, order=3, mode='constant', cval=0.0)
DAPI_1_rotate_shift = shift(DAPI_1_rotate, [-dx, -dy], order=3, mode='constant')

plt.figure(figsize=(16, 10))
plt.subplot(131)
plt.imshow(DAPI_1_processed, cmap='gray')
plt.title('DAPI_1')
plt.axis('off')
plt.subplot(132)
plt.imshow(DAPI_2_processed, cmap='gray')
plt.title('DAPI_2')
plt.axis('off')
plt.subplot(133)
plt.imshow(DAPI_2_processed, cmap='gray')
plt.imshow(DAPI_1_rotate_shift, cmap='twilight', alpha=0.5)
plt.title('Aligned images (DAPI_2 as reference - DAPI_1 in blue)')
plt.axis('off')
plt.show()
