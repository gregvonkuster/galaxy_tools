#!/usr/bin/env python
import argparse
import os
import shutil

from dipy.core.gradients import gradient_table
from dipy.data import fetch_stanford_hardi, fetch_stanford_t1, read_stanford_labels
from dipy.io import read_bvals_bvecs
from dipy.io.image import save_nifti
from matplotlib import pyplot

import nibabel

parser = argparse.ArgumentParser()
parser.add_argument('--drmi_dataset', dest='drmi_dataset', help='Input dataset')
parser.add_argument('--drmi_dataset_type', dest='drmi_dataset_type', help='Input dataset type')
parser.add_argument('--output_nifti1', dest='output_nifti1', help='Output Nifti1 dataset')
parser.add_argument('--output_nifti1_extra_files', dest='output_nifti1_extra_files', help='Output Nifti1 extra files')
parser.add_argument('--output_png', dest='output_png', default=None, help='Output dataset')

args = parser.parse_args()

def move_directory_files(source_dir, destination_dir, copy=False, remove_source_dir=False):
    source_directory = os.path.abspath(source_dir)
    destination_directory = os.path.abspath(destination_dir)
    if not os.path.isdir(destination_directory):
        os.makedirs(destination_directory)
    for dir_entry in os.listdir(source_directory):
        source_entry = os.path.join(source_directory, dir_entry)
        if copy:
            shutil.copy(source_entry, destination_directory)
        else:
            shutil.move(source_entry, destination_directory)
    if remove_source_dir:
        os.rmdir(source_directory)

# Get input data.
input_dir = 'stanford_hardi'
if args.drmi_dataset == 'stanford_hardi':
    if args.drmi_dataset_type == "dataset":
        fetch_stanford_hardi()
        fdwi = os.path.join(input_dir, 'HARDI150.nii.gz')
        fbval = os.path.join(input_dir, 'HARDI150.bval')
        fbvec = os.path.join(input_dir, 'HARDI150.bvec')
        img = nibabel.load(fdwi)
    else:
        img, gtab, labels = read_stanford_labels()
        fdwi = os.path.join(input_dir, 'HARDI150.nii.gz')
        fbval = os.path.join(input_dir, 'HARDI150.bval')
        fbvec = os.path.join(input_dir, 'HARDI150.bvec')
else:
    # args.drmi_dataset == 'stanford_t1
        fetch_stanford_t1()
        fdwi = os.path.join(input_dir, 't1.nii.gz')
        fbval = os.path.join(input_dir, 't1.bval')
        fbvec = os.path.join(input_dir, 't1.bvec')
        img = nibabel.load(fdwi)

data = img.get_data()
# data is a 4D array where the first 3 dimensions are the i, j,
# k voxel coordinates and the last dimension is the number of
# non-weighted (S0s) and diffusion-weighted volumes.

if args.drmi_dataset == 'stanford_hardi':
    # Visualize the results using matplotlib.
    axial_middle = data.shape[2] // 2
    pyplot.subplot(1, 2, 1).set_axis_off()
    pyplot.imshow(data[:, :, axial_middle, 0].T, cmap='gray', origin='lower')
    pyplot.subplot(1, 2, 2).set_axis_off()
    pyplot.imshow(data[:, :, axial_middle, 10].T, cmap='gray', origin='lower')
    pyplot.savefig('data.png', bbox_inches='tight')
    shutil.move('data.png', args.output_png)
    # Load the b-values and b-vectors.
    bvals, bvecs = read_bvals_bvecs(fbval, fbvec)
    gtab = gradient_table(bvals, bvecs)
    # gtab can be used to tell what part of the data is the S0
    # volumes (volumes which correspond to b-values of 0).
    S0s = data[:, :, :, gtab.b0s_mask]
    # Save this in a new Nifti file.
    nibabel.save(nibabel.Nifti1Image(S0s, img.affine), 'output.nii')
else:
    save_nifti('output.nii', data, img.affine)

shutil.move('output.nii', args.output_nifti1)
# Move the entire contents of input_dir to output_nifti1_extra_files.
move_directory_files(input_dir, args.output_nifti1_extra_files)
