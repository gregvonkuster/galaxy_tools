#!/usr/bin/env python
import argparse
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
parser.add_argument('--input', dest='input', help='Track Visualization Header dataset')
parser.add_argument('--output_life_candidates', dest='output_life_candidates', help='Output life candidates')

args = parser.parse_args()

# We'll need to know where the corpus callosum is from these variables.
hardi_img, gtab, labels_img = read_stanford_labels()
labels = labels_img.get_data()
cc_slice = labels == 2
fetch_stanford_t1()
t1 = read_stanford_t1()
t1_data = t1.get_data()
data = hardi_img.get_data()

# Read the candidates from file in voxel space:
candidate_sl = [s[0] for s in nib.trackvis.read(args.input, points_space='voxel')[0]]
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
