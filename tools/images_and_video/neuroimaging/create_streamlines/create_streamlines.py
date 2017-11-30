#!/usr/bin/env python
import argparse
import os
import shutil

from dipy.data import read_stanford_hardi
from dipy.reconst import peaks, shm
from dipy.tracking import utils
from dipy.tracking.eudx import EuDX
from dipy.viz import fvtk
from dipy.viz.colormap import line_colors

import matplotlib.pyplot as plt

import nibabel as nib

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--input', dest='input', help='Input dataset')
parser.add_argument('--input_extra_files_path', dest='input_extra_files_path', help='Input dataset extra files paths')
parser.add_argument('--output_superiorfrontal_nifti', dest='output_superiorfrontal_nifti', help='Output superiorfrontal nifti1 dataset')
parser.add_argument('--output_superiorfrontal_nifti_files_path', dest='output_superiorfrontal_nifti_files_path', help='Output superiorfrontal nifti1 extra files path')
parser.add_argument('--output_trackvis_header', dest='output_trackvis_header', help='Output superiorfrontal track visualization header dataset')

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
# TODO: do not hard-code 'stanford_hardi'
input_dir = 'stanford_hardi'
os.mkdir(input_dir)
for f in os.listdir(args.input_extra_files_path):
    shutil.copy(os.path.join(args.input_extra_files_path, f), input_dir)
hardi_img, gtab = read_stanford_hardi()
data = hardi_img.get_data()
labels_file = os.path.join(input_dir, "aparc-reduced.nii.gz")
labels_img = nib.load(labels_file)
labels = labels_img.get_data()

white_matter = (labels == 1) | (labels == 2)
csamodel = shm.CsaOdfModel(gtab, 6)
csapeaks = peaks.peaks_from_model(model=csamodel, data=data, sphere=peaks.default_sphere, relative_peak_threshold=.8, min_separation_angle=45, mask=white_matter)
seeds = utils.seeds_from_mask(white_matter, density=2)
streamline_generator = EuDX(csapeaks.peak_values, csapeaks.peak_indices, odf_vertices=peaks.default_sphere.vertices, a_low=.05, step_sz=.5, seeds=seeds)
affine = streamline_generator.affine
streamlines = list(streamline_generator)
cc_slice = labels == 2
cc_streamlines = utils.target(streamlines, cc_slice, affine=affine)
cc_streamlines = list(cc_streamlines)

M, grouping = utils.connectivity_matrix(cc_streamlines, labels, affine=affine, return_mapping=True, mapping_as_streamlines=True)
lr_superiorfrontal_track = grouping[11, 54]
shape = labels.shape
dm = utils.density_map(lr_superiorfrontal_track, shape, affine=affine)
# Save density map
dm_img = nib.Nifti1Image(dm.astype("int16"), hardi_img.affine)
dm_img.to_filename("lr-superiorfrontal-dm.nii")
shutil.move('lr-superiorfrontal-dm.nii', args.output_superiorfrontal_nifti)
move_directory_files(input_dir, args.output_superiorfrontal_nifti_files_path)

# Make a trackvis header so we can save streamlines
voxel_size = labels_img.header.get_zooms()
trackvis_header = nib.trackvis.empty_header()
trackvis_header['voxel_size'] = voxel_size
trackvis_header['dim'] = shape
trackvis_header['voxel_order'] = "RAS"
# Move streamlines to "trackvis space"
trackvis_point_space = utils.affine_for_trackvis(voxel_size)
lr_sf_trk = utils.move_streamlines(lr_superiorfrontal_track, trackvis_point_space, input_space=affine)
lr_sf_trk = list(lr_sf_trk)
# Save streamlines
for_save = [(sl, None, None) for sl in lr_sf_trk]
nib.trackvis.write("lr-superiorfrontal.trk", for_save, trackvis_header)
shutil.move('lr-superiorfrontal.trk', args.output_trackvis_header)
