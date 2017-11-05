#!/usr/bin/env python
import argparse
import shutil

from dipy.data import fetch_sherbrooke_3shell
from dipy.data import fetch_stanford_hardi
from dipy.data import get_sphere
from dipy.data import read_sherbrooke_3shell
from dipy.data import read_stanford_hardi
from dipy.reconst import dti
from dipy.reconst.dti import color_fa
from dipy.reconst.dti import fractional_anisotropy
from dipy.segment.mask import median_otsu
from dipy.viz import fvtk
from matplotlib import pyplot

import nibabel
import numpy

# http://nipy.org/dipy/examples_built/reconst_dti.html#example-reconst-dti
parser = argparse.ArgumentParser()
parser.add_argument('--drmi_dataset', dest='drmi_dataset', help='Input dataset')
parser.add_argument('--output_nifti1_fa', dest='output_nifti1_fa', help='Output fractional anisotropy Nifti1 dataset')
parser.add_argument('--output_nifti1_evecs', dest='output_nifti1_evecs', help='Output eigen vectors Nifti1 dataset')
parser.add_argument('--output_nifti1_md', dest='output_nifti1_md', help='Output mean diffusivity Nifti1 dataset')
parser.add_argument('--output_nifti1_rgb', dest='output_nifti1_rgb', help='Output RGB-map Nifti1 dataset')
parser.add_argument('--output_png_ellipsoids', dest='output_png_ellipsoids', help='Output ellipsoids PNG dataset')
parser.add_argument('--output_png_odfs', dest='output_png_odfs', help='Output orientation distribution functions PNG dataset')
parser.add_argument('--output_png_middle_axial_slice', dest='output_png_middle_axial_slice', help='Output middle axial slice PNG dataset')

args = parser.parse_args()

# Get input data.
input_dir = args.drmi_dataset
if input_dir == 'sherbrooke_3shell':
    fetch_sherbrooke_3shell()
    img, gtab = read_sherbrooke_3shell()
elif input_dir == 'stanford_hardi':
    fetch_stanford_hardi()
    img, gtab = read_stanford_hardi()

data = img.get_data()
maskdata, mask = median_otsu(data, 3, 1, True, vol_idx=range(10, 50), dilate=2)

axial_middle = data.shape[2] // 2
pyplot.subplot(1, 2, 1).set_axis_off()
pyplot.imshow(data[:, :, axial_middle, 0].T, cmap='gray', origin='lower')
pyplot.subplot(1, 2, 2).set_axis_off()
pyplot.imshow(data[:, :, axial_middle, 10].T, cmap='gray', origin='lower')
pyplot.savefig('middle_axial.png', bbox_inches='tight')
shutil.move('middle_axial.png', args.output_png_middle_axial_slice)

tenmodel = dti.TensorModel(gtab)
tenfit = tenmodel.fit(maskdata)

fa = fractional_anisotropy(tenfit.evals)
fa[numpy.isnan(fa)] = 0
fa_img = nibabel.Nifti1Image(fa.astype(numpy.float32), img.affine)
nibabel.save(fa_img, 'output_fa.nii')
shutil.move('output_fa.nii', args.output_nifti1_fa)

evecs_img = nibabel.Nifti1Image(tenfit.evecs.astype(numpy.float32), img.affine)
nibabel.save(evecs_img, 'output_evecs.nii')
shutil.move('output_evecs.nii', args.output_nifti1_evecs)

md1 = dti.mean_diffusivity(tenfit.evals)
nibabel.save(nibabel.Nifti1Image(md1.astype(numpy.float32), img.affine), 'output_md.nii')
shutil.move('output_md.nii', args.output_nifti1_md)

fa = numpy.clip(fa, 0, 1)
rgb = color_fa(fa, tenfit.evecs)
nibabel.save(nibabel.Nifti1Image(numpy.array(255 * rgb, 'uint8'), img.affine), 'output_rgb.nii')
shutil.move('output_rgb.nii', args.output_nifti1_rgb)

sphere = get_sphere('symmetric724')
ren = fvtk.ren()

evals = tenfit.evals[13:43, 44:74, 28:29]
evecs = tenfit.evecs[13:43, 44:74, 28:29]
cfa = rgb[13:43, 44:74, 28:29]
cfa /= cfa.max()
fvtk.add(ren, fvtk.tensor(evals, evecs, cfa, sphere))
fvtk.record(ren, n_frames=1, out_path='tensor_ellipsoids.png', size=(600, 600))
shutil.move('tensor_ellipsoids.png', args.output_png_ellipsoids)

fvtk.clear(ren)

tensor_odfs = tenmodel.fit(data[20:50, 55:85, 38:39]).odf(sphere)
fvtk.add(ren, fvtk.sphere_funcs(tensor_odfs, sphere, colormap=None))
fvtk.record(ren, n_frames=1, out_path='tensor_odfs.png', size=(600, 600))
shutil.move('tensor_odfs.png', args.output_png_odfs)
