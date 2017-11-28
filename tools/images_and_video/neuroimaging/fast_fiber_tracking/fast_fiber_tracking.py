#!/usr/bin/env python
import argparse
import shutil

from dipy.data import fetch_sherbrooke_3shell
from dipy.data import fetch_stanford_hardi
from dipy.data import get_sphere
from dipy.data import read_sherbrooke_3shell
from dipy.data import read_stanford_hardi
from dipy.direction import peaks_from_model
from dipy.io.image import save_nifti
from dipy.reconst.csdeconv import (ConstrainedSphericalDeconvModel, auto_response)
from dipy.reconst.dti import TensorModel
from dipy.segment.mask import median_otsu
from dipy.tracking.local import LocalTracking, ThresholdTissueClassifier
from dipy.tracking.streamline import Streamlines
from dipy.tracking.utils import random_seeds_from_mask
from dipy.viz import actor, window

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--drmi_dataset', dest='drmi_dataset', help='Input dataset')
parser.add_argument('--output_csd_direction_field', dest='output_csd_direction_field', help='Output csd direction field dataset')
parser.add_argument('--output_det_streamlines', dest='output_det_streamlines', help='Output det streamlines dataset')
parser.add_argument('--output_fa_map', dest='output_fa_map', help='Output fa map dataset')

args = parser.parse_args()

interactive = False

# Get input data.
input_dir = args.drmi_dataset
if input_dir == 'sherbrooke_3shell':
    fetch_sherbrooke_3shell()
    img, gtab = read_sherbrooke_3shell()
elif input_dir == 'stanford_hardi':
    fetch_stanford_hardi()
    img, gtab = read_stanford_hardi()

data = img.get_data()
maskdata, mask = median_otsu(data, 3, 1, False, vol_idx=range(10, 50), dilate=2)

response, ratio = auto_response(gtab, data, roi_radius=10, fa_thr=0.7)
csd_model = ConstrainedSphericalDeconvModel(gtab, response)
sphere = get_sphere('symmetric724')
csd_peaks = peaks_from_model(model=csd_model, data=data, sphere=sphere, mask=mask, relative_peak_threshold=.5, min_separation_angle=25, parallel=True)

tensor_model = TensorModel(gtab, fit_method='WLS')
tensor_fit = tensor_model.fit(data, mask)
fa = tensor_fit.fa

tissue_classifier = ThresholdTissueClassifier(fa, 0.1)
seeds = random_seeds_from_mask(fa > 0.3, seeds_count=1)

ren = window.Renderer()
ren.add(actor.peak_slicer(csd_peaks.peak_dirs, csd_peaks.peak_values, colors=None))
window.record(ren, out_path='csd_direction_field.png', size=(900, 900))
shutil.move('csd_direction_field.png', args.output_csd_direction_field)

streamline_generator = LocalTracking(csd_peaks, tissue_classifier, seeds, affine=np.eye(4), step_size=0.5)
streamlines = Streamlines(streamline_generator)

ren.clear()
ren.add(actor.line(streamlines))
window.record(ren, out_path='det_streamlines.png', size=(900, 900))
shutil.move('det_streamlines.png', args.output_det_streamlines)

save_nifti('fa_map.nii', fa, img.affine)
shutil.move('fa_map.nii', args.output_fa_map)
