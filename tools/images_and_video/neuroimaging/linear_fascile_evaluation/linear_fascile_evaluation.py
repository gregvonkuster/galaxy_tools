#!/usr/bin/env python
import argparse
import os
import shutil

import dipy.core.optimize as opt
import dipy.tracking.life as life
from dipy.data import fetch_stanford_t1, read_stanford_labels, read_stanford_t1
from dipy.viz import fvtk
from dipy.viz.colormap import line_colors

import matplotlib
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import AxesGrid

import nibabel as nib

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--input_nifti1', dest='input_nifti1', help='Input nifti1 dataset')
parser.add_argument('--input_nifti1_files_path', dest='input_nifti1_files_path', help='Input nifti1 extra files path')
parser.add_argument('--input_nifti2', dest='input_nifti2', help='Input nifti2 dataset')
parser.add_argument('--input_nifti2_files_path', dest='input_nifti2_files_path', help='Input nifti2 extra files path')
parser.add_argument('--input_trackvis', dest='input_trackvis', help='Track Visualization Header dataset')
parser.add_argument('--output_life_candidates', dest='output_life_candidates', help='Output life candidates')

args = parser.parse_args()

# Get input data.
# TODO: do not hard-code 'stanford_hardi'
input_dir = 'stanford_hardi'
os.mkdir(input_dir)
# Copy the dRMI dataset (stanford_t1) files.
for f in os.listdir(args.input_nifti1_files_path):
    shutil.copy(os.path.join(args.input_nifti1_files_path, f), input_dir)
# Copy the dRMI dataset and label map (stanford_hardi) files.
for f in os.listdir(args.input_nifti2_files_path):
    shutil.copy(os.path.join(args.input_nifti2_files_path, f), input_dir)

# We'll need to know where the corpus callosum is from these variables.
hardi_img, gtab, labels_img = read_stanford_labels()
labels = labels_img.get_data()
cc_slice = labels == 2
t1 = read_stanford_t1()
t1_data = t1.get_data()
data = hardi_img.get_data()

# Read the candidates from file in voxel space:
candidate_sl = [s[0] for s in nib.trackvis.read(args.input_trackvis, points_space='voxel')[0]]
# Visualize the initial candidate group of streamlines
# in 3D, relative to the anatomical structure of this brain.
candidate_streamlines_actor = fvtk.streamtube(candidate_sl, line_colors(candidate_sl))
cc_ROI_actor = fvtk.contour(cc_slice, levels=[1], colors=[(1., 1., 0.)], opacities=[1.])
vol_actor = fvtk.slicer(t1_data)
vol_actor.display(40, None, None)
vol_actor2 = vol_actor.copy()
vol_actor2.display(None, None, 35)
# Add display objects to canvas.
ren = fvtk.ren()
fvtk.add(ren, candidate_streamlines_actor)
fvtk.add(ren, cc_ROI_actor)
fvtk.add(ren, vol_actor)
fvtk.add(ren, vol_actor2)
fvtk.record(ren, n_frames=1, out_path="life_candidates.png", size=(800, 800))
shutil.move("life_candidates.png", args.output_life_candidates)
