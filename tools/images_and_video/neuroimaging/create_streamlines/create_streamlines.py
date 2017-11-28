#!/usr/bin/env python
import argparse
import shutil

from dipy.data import fetch_stanford_t1, read_stanford_labels, read_stanford_t1
from dipy.reconst import peaks, shm
from dipy.tracking import utils
from dipy.tracking.eudx import EuDX
from dipy.viz import fvtk
from dipy.viz.colormap import line_colors

import matplotlib.pyplot as plt

import nibabel as nib

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--drmi_dataset', dest='drmi_dataset', help='Input dataset')
parser.add_argument('--output_corpuscallosum_axial', dest='output_corpuscallosum_axial', help='Output corpuscallosum axial dataset')
parser.add_argument('--output_corpuscallosum_sagittal', dest='output_corpuscallosum_sagittal', help='Output corpuscallosum sagittal dataset')
parser.add_argument('--output_connectivity', dest='output_connectivity', help='Output connectivity dataset')
parser.add_argument('--output_superiorfrontal_nifti', dest='output_superiorfrontal_nifti', help='Output superiorfrontal nifti1 dataset')
parser.add_argument('--output_trackvis_header', dest='output_trackvis_header', help='Output superiorfrontal track visualization header dataset')

args = parser.parse_args()

hardi_img, gtab, labels_img = read_stanford_labels()
data = hardi_img.get_data()
labels = labels_img.get_data()

# For possible later use: if args.drmi_dataset == 'stanford_t1':
fetch_stanford_t1()
t1 = read_stanford_t1()
t1_data = t1.get_data()
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
other_streamlines = utils.target(streamlines, cc_slice, affine=affine, include=False)
other_streamlines = list(other_streamlines)
assert len(other_streamlines) + len(cc_streamlines) == len(streamlines)
# Make display objects
color = line_colors(cc_streamlines)
cc_streamlines_actor = fvtk.line(cc_streamlines, line_colors(cc_streamlines))
cc_ROI_actor = fvtk.contour(cc_slice, levels=[1], colors=[(1., 1., 0.)], opacities=[1.])
vol_actor = fvtk.slicer(t1_data)
vol_actor.display(40, None, None)
vol_actor2 = vol_actor.copy()
vol_actor2.display(None, None, 35)
# Add display objects to canvas
r = fvtk.ren()
fvtk.add(r, vol_actor)
fvtk.add(r, vol_actor2)
fvtk.add(r, cc_streamlines_actor)
fvtk.add(r, cc_ROI_actor)
# Save figures
fvtk.record(r, n_frames=1, out_path="corpuscallosum_axial.png", size=(800, 800))
shutil.move("corpuscallosum_axial.png", args.output_corpuscallosum_axial)
fvtk.camera(r, [-1, 0, 0], [0, 0, 0], viewup=[0, 0, 1])
fvtk.record(r, n_frames=1, out_path="corpuscallosum_sagittal.png", size=(800, 800))
shutil.move("corpuscallosum_sagittal.png", args.output_corpuscallosum_sagittal)
M, grouping = utils.connectivity_matrix(cc_streamlines, labels, affine=affine, return_mapping=True, mapping_as_streamlines=True)
M[:3, :] = 0
M[:, :3] = 0
plt.imshow(np.log1p(M), interpolation='nearest')
plt.savefig("connectivity.png")
shutil.move("connectivity.png", args.output_connectivity)
lr_superiorfrontal_track = grouping[11, 54]
shape = labels.shape
dm = utils.density_map(lr_superiorfrontal_track, shape, affine=affine)
# Save density map
dm_img = nib.Nifti1Image(dm.astype("int16"), hardi_img.affine)
dm_img.to_filename("lr-superiorfrontal-dm.nii")
shutil.move('lr-superiorfrontal-dm.nii', args.output_superiorfrontal_nifti)
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
