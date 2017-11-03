#!/usr/bin/env python
import argparse
import os
import nibabel
import shutil

from dipy.core.gradients import gradient_table
from dipy.data import fetch_sherbrooke_3shell
from dipy.io import read_bvals_bvecs
from matplotlib import pyplot

parser = argparse.ArgumentParser()
parser.add_argument('--input', dest='input', help='Input dataset')
parser.add_argument('--output_nifti1', dest='output_nifti1', help='Output Nifti1 dataset')
parser.add_argument('--output_png', dest='output_png', help='Output dataset')

args = parser.parse_args()

input_dir = 'sherbrooke_3shell'
# Get input data.
fetch_sherbrooke_3shell()
fdwi = os.path.join(input_dir, 'HARDI193.nii.gz')
fbval = os.path.join(input_dir, 'HARDI193.bval')
fbvec = os.path.join(input_dir, 'HARDI193.bvec')
# Load the dMRI datasets.
img = nibabel.load(fdwi)
data = img.get_data()
# data is a 4D array where the first 3 dimensions are the i, j,
# k voxel coordinates and the last dimension is the number of
# non-weighted (S0s) and diffusion-weighted volumes.
# Visualize the results using matplotlib.
axial_middle = data.shape[2] // 2
pyplot.figure('Showing the datasets')
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
shutil.move('output.nii', args.output_nifti1)
