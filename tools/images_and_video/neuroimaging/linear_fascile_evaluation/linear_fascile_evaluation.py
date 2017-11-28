#!/usr/bin/env python
import argparse
import numpy as np
import os.path as op
import nibabel as nib
import dipy.core.optimize as opt
import dipy.tracking.life as life
import matplotlib.pyplot as plt
import matplotlib

from dipy.viz.colormap import line_colors
from dipy.viz import fvtk
from mpl_toolkits.axes_grid1 import AxesGrid
from dipy.data import read_stanford_labels, fetch_stanford_t1, read_stanford_t1

parser = argparse.ArgumentParser()
parser.add_argument('--candidates', dest='candidates', help='Candidates selection')
parser.add_argument('--output_life_candidates', dest='output_life_candidates', help='Output life candidates')
parser.add_argument('--output_life_optimized', dest='output_life_optimized', help='Output life optimized streamlines')
parser.add_argument('--output_beta_histogram', dest='output_beta_histogram', help='Output beta histogram')
parser.add_argument('--output_error_histograms', dest='output_error_histograms', help='Output error histograms')
parser.add_argument('--output_spatial_errors', dest='output_spatial_errors', help='Output spatial errors')

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
candidate_sl = [s[0] for s in nib.trackvis.read(args.candidates, points_space='voxel')[0]]
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
# Initialize a LiFE model.
fiber_model = life.FiberModel(gtab)
# Fit the model, producing a FiberFit class instance,
# that stores the data, as well as the results of the
# fitting procedure.
fiber_fit = fiber_model.fit(data, candidate_sl, affine=np.eye(4))
fig, ax = plt.subplots(1)
ax.hist(fiber_fit.beta, bins=100, histtype='step')
ax.set_xlabel('Fiber weights')
ax.set_ylabel('# fibers')
fig.savefig("beta_histogram.png")
shutil.move("beta_histogram.png", args.output_beta_histogram)
# Filter out these redundant streamlines and
# generate an optimized group of streamlines.
optimized_sl = list(np.array(candidate_sl)[np.where(fiber_fit.beta>0)[0]])
ren = fvtk.ren()
fvtk.add(ren, fvtk.streamtube(optimized_sl, line_colors(optimized_sl)))
fvtk.add(ren, cc_ROI_actor)
fvtk.add(ren, vol_actor)
fvtk.record(ren, n_frames=1, out_path="optimized.png", size=(800, 800))
shutil.move("optimized.png", args.output_life_optimized)
model_predict = fiber_fit.predict()
# Focus on the error in prediction of the diffusion-weighted
# data, and calculate the root of the mean squared error.
model_error = model_predict - fiber_fit.data
model_rmse = np.sqrt(np.mean(model_error[:, 10:] ** 2, -1))
# Calculate another error term by assuming that the weight for each streamline
# is equal to zero. This produces the naive prediction of the mean of the
# signal in each voxel.
beta_baseline = np.zeros(fiber_fit.beta.shape[0])
pred_weighted = np.reshape(opt.spdot(fiber_fit.life_matrix, beta_baseline), (fiber_fit.vox_coords.shape[0], np.sum(~gtab.b0s_mask)))
mean_pred = np.empty((fiber_fit.vox_coords.shape[0], gtab.bvals.shape[0]))
S0 = fiber_fit.b0_signal
# Since the fitting is done in the demeaned S/S0 domain,
# add back the mean and then multiply by S0 in every voxel:
mean_pred[..., gtab.b0s_mask] = S0[:, None]
mean_pred[..., ~gtab.b0s_mask] = (pred_weighted + fiber_fit.mean_signal[:, None]) * S0[:, None]
mean_error = mean_pred - fiber_fit.data
mean_rmse = np.sqrt(np.mean(mean_error ** 2, -1))
# Compare the overall distribution of errors
# between these two alternative models of the ROI.
fig, ax = plt.subplots(1)
ax.hist(mean_rmse - model_rmse, bins=100, histtype='step')
ax.text(0.2, 0.9,'Median RMSE, mean model: %.2f' % np.median(mean_rmse), horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
ax.text(0.2, 0.8,'Median RMSE, LiFE: %.2f' % np.median(model_rmse), horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
ax.set_xlabel('RMS Error')
ax.set_ylabel('# voxels')
fig.savefig("error_histograms.png")
shutil.move("error_histograms.png", args.output_error_histograms)
# Show the spatial distribution of the two error terms,
# and of the improvement with the model fit:
vol_model = np.ones(data.shape[:3]) * np.nan
vol_model[fiber_fit.vox_coords[:, 0], fiber_fit.vox_coords[:, 1], fiber_fit.vox_coords[:, 2]] = model_rmse
vol_mean = np.ones(data.shape[:3]) * np.nan
vol_mean[fiber_fit.vox_coords[:, 0], fiber_fit.vox_coords[:, 1], fiber_fit.vox_coords[:, 2]] = mean_rmse
vol_improve = np.ones(data.shape[:3]) * np.nan
vol_improve[fiber_fit.vox_coords[:, 0], fiber_fit.vox_coords[:, 1], fiber_fit.vox_coords[:, 2]] = mean_rmse - model_rmse
sl_idx = 49
fig = plt.figure()
fig.subplots_adjust(left=0.05, right=0.95)
ax = AxesGrid(fig, 111, nrows_ncols = (1, 3), label_mode = "1", share_all = True, cbar_location="top", cbar_mode="each", cbar_size="10%", cbar_pad="5%")
ax[0].matshow(np.rot90(t1_data[sl_idx, :, :]), cmap=matplotlib.cm.bone)
im = ax[0].matshow(np.rot90(vol_model[sl_idx, :, :]), cmap=matplotlib.cm.hot)
ax.cbar_axes[0].colorbar(im)
ax[1].matshow(np.rot90(t1_data[sl_idx, :, :]), cmap=matplotlib.cm.bone)
im = ax[1].matshow(np.rot90(vol_mean[sl_idx, :, :]), cmap=matplotlib.cm.hot)
ax.cbar_axes[1].colorbar(im)
ax[2].matshow(np.rot90(t1_data[sl_idx, :, :]), cmap=matplotlib.cm.bone)
im = ax[2].matshow(np.rot90(vol_improve[sl_idx, :, :]), cmap=matplotlib.cm.RdBu)
ax.cbar_axes[2].colorbar(im)
for lax in ax:
    lax.set_xticks([])
    lax.set_yticks([])
fig.savefig("spatial_errors.png")
shutil.move("spatial_errors.png", args.output_spatial_errors)
