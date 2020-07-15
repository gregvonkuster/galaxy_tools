#!/usr/bin/env python

import numpy as np
from numpy.random import rand, randn

import metrics


def gen_sphere_grid(data, mles_param_dict, prior, param_dict, scaleFactor):
    x = mles_param_dict['xc']
    y = mles_param_dict['yc']
    z = mles_param_dict['zc']

    if rand() < 0.2:
        px = x + param_dict['xTau'] * randn()
        py = y + param_dict['yTau'] * randn()
        pz = z + param_dict['zTau'] * randn()
    else:
        px,  py,  pz = x,  y,  z

    xE = mles_param_dict['xE']
    yE = mles_param_dict['yE']
    zE = mles_param_dict['zE']
    r = mles_param_dict['r']
    if rand() < 0.5:
        if rand() < 0.7:
            scale = abs(mles_param_dict['rSig'] * randn())/1000
        elif rand() < 0.7:
            scale = abs(mles_param_dict['rSig'] * randn())/100
        else:
            scale = abs(mles_param_dict['rSig'] * randn())/10

        scale = scale * scaleFactor

        gr = abs(r + r * (1-mles_param_dict['fix_r']) * scale * param_dict['grid'])
        gxE = abs(xE + xE * (1-mles_param_dict['fix_xE']) * scale * param_dict['grid'])
        gyE = abs(yE + yE * (1-mles_param_dict['fix_yE']) * scale * param_dict['grid'])
        gzE = abs(zE + zE * (1-mles_param_dict['fix_zE']) * scale * param_dict['grid'])
    else:
        if rand() < 0.7:
            scale = abs(mles_param_dict['rSig'] * randn())/1000
        elif rand() < 0.7:
            scale = abs(mles_param_dict['rSig'] * randn())/100
        else:
            scale = abs(mles_param_dict['rSig'] * randn())/10

        scale = scale * scaleFactor
        gr = abs(r + r * (1-mles_param_dict['fix_r']) * scale * param_dict['grid'])
        if rand() < 0.7:
            scale = abs(mles_param_dict['rSig'] * randn())/1000
        elif rand() < 0.7:
            scale = abs(mles_param_dict['rSig'] * randn())/100
        else:
            scale = abs(mles_param_dict['rSig'] * randn())/10

        scale = scale * scaleFactor
        gxE = abs(xE + xE * (1-mles_param_dict['fix_xE']) * scale * param_dict['grid'])
        if rand() < 0.7:
            scale = abs(mles_param_dict['rSig'] * randn())/1000
        elif rand() < 0.7:
            scale = abs(mles_param_dict['rSig'] * randn())/100
        else:
            scale = abs(mles_param_dict['rSig'] * randn())/10

        scale = scale * scaleFactor
        gyE = abs(yE + yE * (1-mles_param_dict['fix_yE']) * scale * param_dict['grid'])
        if rand() < 0.7:
            scale = abs(mles_param_dict['rSig'] * randn())/1000
        elif rand() < 0.7:
            scale = abs(mles_param_dict['rSig'] * randn())/100
        else:
            scale = abs(mles_param_dict['rSig'] * randn())/10

        scale = scale * scaleFactor
        gzE = abs(zE + zE * (1-mles_param_dict['fix_zE']) * scale * param_dict['grid'])

    gLogLike = np.zeros(param_dict['nGrid']**4, )
    for ir in range(param_dict['nGrid']):
        for ixE in range(param_dict['nGrid']):
            for iyE in range(param_dict['nGrid']):
                for izE in range(param_dict['nGrid']):
                    # tMaxIndex = (ir-1) * param_dict['nGrid']^3+(ixE-1) * param_dict['nGrid']^2+(iyE-1) * param_dict['nGrid'] +izE
                    # [izE iyE ixE ir tMaxIndex]
                    gLogLike[(ir) * param_dict['nGrid']**3 + (ixE) * param_dict['nGrid']**2+(iyE) * param_dict['nGrid'] + izE] = metrics.calc_log_like_sphere_mix(px, py, pz, gxE[ixE], gyE[iyE], gzE[izE], gr[ir], data[mles_param_dict['sInGIndex'], :], mles_param_dict['rSig'])
                    # [gr(ir) gxE(ixE) gyE(iyE) gzE(izE) gLogLike(tMaxIndex)]

    maxIndex = np.where(gLogLike == np.amax(gLogLike))
    maxIndex = maxIndex[0][0]
    maxIr = int(np.ceil(maxIndex/(param_dict['nGrid']**3)))
    maxIxE = int(np.ceil((maxIndex-(maxIr-1) * param_dict['nGrid']**3)/(param_dict['nGrid']**2)))
    maxIyE = int(np.ceil((maxIndex-(maxIr-1) * param_dict['nGrid']**3-(maxIxE-1) * param_dict['nGrid']**2)/(param_dict['nGrid'])))
    maxIzE = maxIndex - (maxIr - 1) * param_dict['nGrid']**3 - (maxIxE - 1) * param_dict['nGrid']**2 - (maxIyE - 1) * param_dict['nGrid']
    # tMaxIndex = (maxIr - 1) * param_dict['nGrid']^3+(maxIxE - 1) * param_dict['nGrid']^2+(maxIyE - 1) * param_dict['nGrid'] +maxIzE
    pr = gr[maxIr - 1]
    pxE = gxE[maxIxE - 1]
    pyE = gyE[maxIyE - 1]
    pzE = gzE[maxIzE - 1]
    # tlogLike = metrics.calc_log_like_sphere_mix(px, py, pz, pxE, pyE, pzE, pr, data(mles_param_dict.sInGIndex, :), mles_param_dict['rSig'])

    logLike = metrics.calc_log_like_sphere_mix(x, y, z, xE, yE, zE, r, data[mles_param_dict['sInGIndex'], :], mles_param_dict['rSig'])
    pLogLike = gLogLike[maxIndex]
    logPrior = metrics.calc_log_prior_sphere(x, y, z, r, xE, yE, zE, prior)
    pLogPrior = metrics.calc_log_prior_sphere(px, py, pz, pr, pxE, pyE, pzE, prior)

    lap = pLogLike + pLogPrior - logLike - logPrior

    if np.log(rand()) < lap:
        x = px
        y = py
        z = pz
        r = pr
        xE = pxE
        yE = pyE
        zE = pzE

    return (x, y, z, r, xE, yE, zE)
