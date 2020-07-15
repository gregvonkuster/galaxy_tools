#!/usr/bin/env python

import gen_r_sig_3d
import gen_selection_in_g_3d
import gen_sphere
import gen_sphere_grid
import numpy as np
import matplotlib.pyplot as plt

from copy import copy
from data_utils import map_to_non_euclid_3d
from metrics import calc_log_like_sphere_mix
from mle_priors_3d import MLEPriors3D
from numpy.random import rand
from plot_3d import plot_3d

"""
Maximum Liklihood Estimation for a sphere to fit the data.
In the comments below, 'data points' may be referred as DP.
Spherical Coordinates may be referred as SC.
Cartesian Coordinates may be referred as CC.
"""


def mle_sphere(data, cell_id, param_dict, log_fh):
    n = data.shape[0]

    # Param initialization.
    mles_param_dict = {}
    mles_param_dict['xc'] = np.mean(data[:, 0])
    mles_param_dict['yc'] = np.mean(data[:, 1])
    mles_param_dict['zc'] = np.mean(data[:, 2])
    mles_param_dict['r'] = np.mean(np.sqrt((data[:, 0] - mles_param_dict['xc'])**2 + (data[:, 1] - mles_param_dict['yc'])**2 + (data[:, 2] - mles_param_dict['zc'])**2))
    mles_param_dict['xE'] = 1
    mles_param_dict['yE'] = 1
    mles_param_dict['zE'] = 1
    log_fh.write("Param Init {}\n".format(mles_param_dict))

    # Convert the data to spherical coordinate system. rTheta
    # has 3 axes: Theta (azimuthal angle), Phi (polar angle)
    # and r (radius).
    rTheta = map_to_non_euclid_3d(data,
                                  mles_param_dict['xc'],
                                  mles_param_dict['yc'],
                                  mles_param_dict['zc'],
                                  mles_param_dict['xE'],
                                  mles_param_dict['yE'],
                                  mles_param_dict['zE'])

    # mean radius of DP in spherical coordinates
    meanR = np.mean(rTheta[:, 2])

    # Standard Deviation of radius of DP in spherical coordinates.
    stdR = np.std(rTheta[:, 2])

    # Normalization of radius of DP in spherical coordinates.
    zR = (rTheta[:, 2]-meanR)/stdR
    log_fh.write("meanR {}, stdR {}\n".format(meanR, stdR))

    # Inner points: index of all the normalized data pointsthat are
    # within 1 unit of standard deviation. (This should account to
    # 68.2% of DP)
    mles_param_dict['sInGIndex'] = np.where(abs(zR) < 1)[0]

    # Outer points: index of all the normalized data points that are outside 1 unit of standard deviation.
    mles_param_dict['sOutGIndex'] = np.where(abs(zR) >= 1)[0]

    # log_fh.write("sInGIndex {}, sOutGIndex {}\n".format(mles_param_dict['sInGIndex'], mles_param_dict['sOutGIndex']))
    mles_param_dict['sInG'] = np.zeros((n, ))

    # Marks all the inner points as 1. All the points lying outside are marked as 0.
    mles_param_dict['sInG'][mles_param_dict['sInGIndex']] = 1
    # log_fh.write("sInG {}\n".format(mles_param_dict['sInG']))

    if param_dict['fixR']:
        mles_param_dict['r'] = 1
    else:
        mles_param_dict['r'] = np.mean(rTheta[mles_param_dict['sInGIndex'], 2])

    # Standard Deviation of DP in SC having normalized radius within
    # 1 standard deviation.
    mles_param_dict['rSig'] = np.std(rTheta[mles_param_dict['sInGIndex'], 2])

    # mean of x-coord of DP in CC having normalized radius within 1 standard deviation.
    mles_param_dict['xc'] = np.mean(data[mles_param_dict['sInGIndex'], 0])

    # mean of y-coord of DP in CC having normalized radius within 1 standard deviation.
    mles_param_dict['yc'] = np.mean(data[mles_param_dict['sInGIndex'], 1])

    # mean of z-coord of DP in CC having normalized radius within 1 standard deviation.
    mles_param_dict['zc'] = np.mean(data[mles_param_dict['sInGIndex'], 2])

    sX = np.std(data[mles_param_dict['sInGIndex'], 0])
    sY = np.std(data[mles_param_dict['sInGIndex'], 1])
    sZ = np.std(data[mles_param_dict['sInGIndex'], 2])
    mS = max(sX, max(sY, sZ))

    rTheta = map_to_non_euclid_3d(data,
                                  mles_param_dict['xc'],
                                  mles_param_dict['yc'],
                                  mles_param_dict['zc'],
                                  mles_param_dict['xE'],
                                  mles_param_dict['yE'],
                                  mles_param_dict['zE'])

    tprior = MLEPriors3D(cMean=[mles_param_dict['xc'],
                         mles_param_dict['yc'],
                         mles_param_dict['zc']],
                         cStd=[0.1, 0.1, 0.1],
                         rMean=mles_param_dict['r'],
                         rStd=0.1,
                         eMean=[1, 1, 1],
                         eStd=[0.5, 0.5, 0.5])

    mles_param_dict['fix_xE'] = 0
    mles_param_dict['fix_yE'] = 0
    mles_param_dict['fix_zE'] = 0
    mles_param_dict['fix_r'] = 1

    if not param_dict['fixR']:
        mles_param_dict['fix_r'] = 0
        if mS == sX:
            mles_param_dict['fix_xE'] = 1
        if mS == sY:
            mles_param_dict['fix_yE'] = 1
        if mS == sZ:
            mles_param_dict['fix_zE'] = 1

    mles_param_dict['xE'] = sX/mS
    mles_param_dict['yE'] = sY/mS
    mles_param_dict['zE'] = sZ/mS

    tmp_param_dict = copy(param_dict)
    tmp_param_dict['xETau'] = 10 * param_dict['xETau']
    tmp_param_dict['yETau'] = 10 * param_dict['yETau']
    tmp_param_dict['zETau'] = 10 * param_dict['zETau']
    tmp_param_dict['xTau'] = 0.1 * param_dict['xTau']
    tmp_param_dict['yTau'] = 0.1 * param_dict['yTau']
    tmp_param_dict['zTau'] = 0.1 * param_dict['zTau']

    rTheta = map_to_non_euclid_3d(data,
                                  mles_param_dict['xc'],
                                  mles_param_dict['yc'],
                                  mles_param_dict['zc'],
                                  mles_param_dict['xE'],
                                  mles_param_dict['yE'],
                                  mles_param_dict['zE'])

    tprior = MLEPriors3D(cMean=[mles_param_dict['xc'],
                         mles_param_dict['yc'],
                         mles_param_dict['zc']],
                         cStd=[0.1, 0.1, 0.1],
                         rMean=mles_param_dict['r'],
                         rStd=0.1,
                         eMean=[1, 1, 1],
                         eStd=[0.5, 0.5, 0.5])

    if param_dict['fixR']:
        mles_param_dict['r'] = 1
    else:
        mles_param_dict['r'] = np.mean(rTheta[mles_param_dict['sInGIndex'], 2])

    mles_param_dict['rSig'] = np.std(rTheta[mles_param_dict['sInGIndex'], 2])
    mles_param_dict['xc'] = np.mean(data[mles_param_dict['sInGIndex'], 0])
    mles_param_dict['yc'] = np.mean(data[mles_param_dict['sInGIndex'], 1])
    mles_param_dict['zc'] = np.mean(data[mles_param_dict['sInGIndex'], 2])

    scaleFactor = 10

    for n in range(param_dict['burnin'] // 2):
        log_fh.write("MLE Sphere Loop {}\n".format(n+1))
        tup = gen_sphere.gen_sphere(data, mles_param_dict, tprior, tmp_param_dict)
        mles_param_dict['xc'], mles_param_dict['yc'], mles_param_dict['zc'], mles_param_dict['r'], mles_param_dict['xE'], mles_param_dict['yE'], mles_param_dict['zE'] = tup

        tup = gen_sphere_grid.gen_sphere_grid(data, mles_param_dict, tprior, tmp_param_dict, scaleFactor)
        mles_param_dict['xc'], mles_param_dict['yc'], mles_param_dict['zc'], mles_param_dict['r'], mles_param_dict['xE'], mles_param_dict['yE'], mles_param_dict['zE'] = tup

        if rand() < 0.2:
            mles_param_dict['rSig'] = gen_r_sig_3d.gen_r_sig_3d(data, mles_param_dict, tprior)

            tup = gen_selection_in_g_3d.gen_selection_in_g_3d(data, mles_param_dict, tprior)
            mles_param_dict['sInG'], mles_param_dict['sInGIndex'], mles_param_dict['sOutGIndex'] = tup

            if param_dict['plotFigures']:
                plt.close('all')
                plot_3d(tmp_param_dict, data, mles_param_dict)
                mles_param_dict['cLogLike'] = calc_log_like_sphere_mix(mles_param_dict['xc'],
                                                                       mles_param_dict['yc'],
                                                                       mles_param_dict['zc'],
                                                                       mles_param_dict['xE'],
                                                                       mles_param_dict['yE'],
                                                                       mles_param_dict['zE'],
                                                                       mles_param_dict['r'],
                                                                       data[mles_param_dict['sInGIndex'], :],
                                                                       mles_param_dict['rSig'])
                log_fh.write("cLogLike {}\n".format(mles_param_dict['cLogLike']))
                log_fh.write("Sphere Center and axes length {}\n".format([mles_param_dict['xc'], mles_param_dict['yc'], mles_param_dict['zc'], mles_param_dict['r'], mles_param_dict['xE'], mles_param_dict['yE'], mles_param_dict['zE']]))

    mles_param_dict['cLogLike'] = calc_log_like_sphere_mix(mles_param_dict['xc'],
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

    if param_dict['plotFigures']:
        plot_3d(tmp_param_dict, data, mles_param_dict)

    return mles_param_dict
