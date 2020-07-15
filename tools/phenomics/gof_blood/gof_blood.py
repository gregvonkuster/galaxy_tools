#!/usr/bin/env python

import argparse
import numpy as np
import os

import data_utils
import mle_sphere
import gen_sphere
import gen_sphere_grid
import gen_r_sig_3d
import gen_selection_in_g_3d
import metrics
import param_ss
import mle_priors_3d

DEFAULT_MLE_SPHERE_PARAM_DICT = dict(xc=0, yc=0, zc=0, r=1, rSig=0.3, xE=1.2, yE=1, zE=0.8)
GRID = np.array([-3, -2, -1, 1, 2, 3])
SCALEFACTOR = 1


def get_bounding_box_dict():
    bounding_box_dict = {}
    # xMin, xMax
    bounding_box_dict['x'] = [-2, 2]
    # yMin, yMax
    bounding_box_dict['y'] = [-2, 2]
    # zMin, zMax
    bounding_box_dict['z'] = [-2, 2]
    return bounding_box_dict


def get_base_file_name(file_path):
    base_file_name = os.path.basename(file_path)
    if base_file_name.find(".") > 0:
        # Eliminate the extension.
        return os.path.splitext(base_file_name)[0]
    else:
        return base_file_name


def get_param_dict(args):
    param_dict = dict(burnin=args.burnin,
                      fixR=args.fixR,
                      grid=GRID,
                      nGrid=args.nGrid,
                      plotFigures=args.plotFigures,
                      randomStart=args.randomStart,
                      rTau=args.rTau,
                      sample=args.sample,
                      thin=args.thin,
                      xETau=args.xETau,
                      xTau=args.xTau,
                      yETau=args.yETau,
                      yTau=args.yTau,
                      zETau=args.zETau,
                      zTau=args.zTau)
    return param_dict


parser = argparse.ArgumentParser()

parser.add_argument('--burnin', action='store', dest='burnin', type=int, required=False, default=250, help='Not sure')
parser.add_argument('--input_dir', action='store', dest='input_dir', help='Directory of input csv image point files')
parser.add_argument('--nGrid', action='store', dest='nGrid', type=int, required=False, default=6, help='Not sure')
parser.add_argument('--output_csv_dir', action='store', dest='output_csv_dir', help='Directory to output csv files')
parser.add_argument('--output_json_dir', action='store', dest='output_json_dir', help='Directory to output json files')
parser.add_argument('--output_log', action='store', dest='output_log', help='Output process log file')
parser.add_argument('--plotFigures', action='store', dest='plotFigures', type=bool, required=False, default=False, help='Plot figures')
parser.add_argument('--randomStart', action='store', dest='randomStart', type=bool, required=False, default=True, help='Do not use DEFAULT_MLE_SPHERE_PARAM_DICT')
parser.add_argument("--fixR", action='store', dest='fixR', type=bool, required=False, default=True, help="Not sure")
parser.add_argument('--sample', action='store', dest='sample', type=int, required=False, default=250, help='Not sure')
parser.add_argument('--thin', action='store', dest='thin', type=int, required=False, default=10, help='Not sure')
parser.add_argument('--xTau', action='store', dest='xTau', type=float, required=False, default=0.01, help='Not sure')
parser.add_argument('--yTau', action='store', dest='yTau', type=float, required=False, default=0.01, help='Not sure')
parser.add_argument('--zTau', action='store', dest='zTau', type=float, required=False, default=0.01, help='Not sure')
parser.add_argument('--rTau', action='store', dest='rTau', type=float, required=False, default=0.01, help='Not sure')
parser.add_argument('--xETau', action='store', dest='xETau', type=float, required=False, default=0.01, help='Not sure')
parser.add_argument('--yETau', action='store', dest='yETau', type=float, required=False, default=0.01, help='Not sure')
parser.add_argument('--zETau', action='store', dest='zETau', type=float, required=False, default=0.01, help='Not sure')

args = parser.parse_args()
param_dict = get_param_dict(args)

# TODO: nLevelSet = 150

log_fh = open(args.output_log, "w")

for file_name in sorted(os.listdir(args.input_dir)):
    file_path = os.path.abspath(os.path.join(args.input_dir, file_name))
    # Extract the cell_id from the file name.
    cell_id = get_base_file_name(file_name)

    # Load or generate synthetic level set data
    # if params_dict['genData']:
    if False:
        # TODO: data = genSphereLevelSet(DEFAULT_MLE_SPHERE_PARAM_DICT, bounding_box, param_dict, nLevelSet)
        meanR = [0, 0, 0]
        # TODO: save(PrjCtrl.inputFileLevelSetDataM,'data')
    else:
        data = data_utils.load_raw3d_data(file_path)
        log_fh.write("First 5 data points before normalization: {}\n".format(data[:5]))
        data, meanR = data_utils.normalize_data(data, log_fh)
        log_fh.write("\nFirst 5 data points after normalization: {}\n\n".format(data[:5]))
        log_fh.write("mean radius {}".format(meanR))
        # TODO: nLevelSet = data.shape[0]

    # Set summary parameters
    param_ss = param_ss.ParamSS(data.shape[0], meanR)

    # Starting value for parameters
    if param_dict['randomStart']:
        mles_param_dict = mle_sphere.mle_sphere(data, cell_id, param_dict, log_fh)
    else:
        mles_param_dict = DEFAULT_MLE_SPHERE_PARAM_DICT

    # Set Priors
    prior = mle_priors_3d.MLEPriors3D(cMean=[mles_param_dict['xc'],
                                      mles_param_dict['yc'],
                                      mles_param_dict['zc']],
                                      cStd=[1, 1, 1],
                                      rMean=mles_param_dict['r'],
                                      rStd=0.1,
                                      eMean=[mles_param_dict['xE'],
                                      mles_param_dict['yE'],
                                      mles_param_dict['zE']],
                                      eStd=[1, 1, 1])

    # MCMC Analysis
    for n in range(args.burnin+args.sample):
        if (np.mod(n, args.thin) == 0 or n == 0):
            log_fh.write("\nn {}\n".format(n))
            mles_param_dict['cLogLike'] = metrics.calc_log_like_sphere_mix(mles_param_dict['xc'],
                                                                           mles_param_dict['yc'],
                                                                           mles_param_dict['zc'],
                                                                           mles_param_dict['xE'],
                                                                           mles_param_dict['yE'],
                                                                           mles_param_dict['zE'],
                                                                           mles_param_dict['r'],
                                                                           data[mles_param_dict['sInGIndex'], :],
                                                                           mles_param_dict['rSig'])
            log_fh.write("\nmles_param_dict:\n")
            log_fh.write("cLogLike {}\n".format(mles_param_dict['cLogLike']))
            log_fh.write("xc {}\n".format(mles_param_dict['xc']))
            log_fh.write("yc {}\n".format(mles_param_dict['yc']))
            log_fh.write("zc {}\n".format(mles_param_dict['zc']))
            log_fh.write("r {}\n".format(mles_param_dict['r']))
            log_fh.write("xE {}\n".format(mles_param_dict['xE']))
            log_fh.write("yE {}\n".format(mles_param_dict['yE']))
            log_fh.write("zE {}\n".format(mles_param_dict['zE']))
            log_fh.write("rSig {}\n".format(mles_param_dict['rSig']))

            if args.plotFigures:
                pass

        tup = gen_sphere.gen_sphere(data, mles_param_dict, prior, param_dict)
        mles_param_dict['xc'], mles_param_dict['yc'], mles_param_dict['zc'], mles_param_dict['r'], mles_param_dict['xE'], mles_param_dict['yE'], mles_param_dict['zE'] = tup
        tup = gen_sphere_grid.gen_sphere_grid(data, mles_param_dict, prior, param_dict, SCALEFACTOR)
        mles_param_dict['xc'], mles_param_dict['yc'], mles_param_dict['zc'], mles_param_dict['r'], mles_param_dict['xE'], mles_param_dict['yE'], mles_param_dict['zE'] = tup
        mles_param_dict['rSig'] = gen_r_sig_3d.gen_r_sig_3d(data, mles_param_dict, prior)
        mles_param_dict['sInG'], mles_param_dict['sInGIndex'], mles_param_dict['sOutGIndex'] = gen_selection_in_g_3d.gen_selection_in_g_3d(data, mles_param_dict, prior)

        if n > args.burnin:
            param_ss.set_params(mles_param_dict)
            # param_ss = storeParam3D(ParamSS, param_dict)

    log_fh.close()

    # Summarize Parameter and print reports
    param_ss.summarize_params(args.sample)
    param_ss.output_csv(args.output_csv_dir, cell_id)
    param_ss.output_json(args.output_json_dir, cell_id)
