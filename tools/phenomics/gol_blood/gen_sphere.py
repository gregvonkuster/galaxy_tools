#!/usr/bin/env python

import numpy as np
from numpy.random import rand,  randn
import metrics


def gen_sphere(data,  mles_param_dict,  prior, param_dict):
    x = mles_param_dict['xc']
    y = mles_param_dict['yc']
    z = mles_param_dict['zc']
    xE = mles_param_dict['xE']
    yE = mles_param_dict['yE']
    zE = mles_param_dict['zE']
    r = mles_param_dict['r']

    # preturb circle,  x
    px = x + param_dict['xTau'] * randn()
    py = y
    pz = z
    pr = r
    pxE = xE
    pyE = yE
    pzE = zE

    logLike = metrics.calc_log_like_sphere_mix(x, y, z, xE, yE, zE, r,  data[mles_param_dict['sInGIndex'], :],  mles_param_dict['rSig'])
    pLogLike = metrics.calc_log_like_sphere_mix(px, py, pz, pxE, pyE, pzE, pr,  data[mles_param_dict['sInGIndex'], :],  mles_param_dict['rSig'])
    logPrior = metrics.calc_log_prior_sphere(x, y, z, r, xE, yE, zE, prior)
    pLogPrior = metrics.calc_log_prior_sphere(px, py, pz, pxE, pyE, pzE, pr, prior)

    lap = pLogLike + pLogPrior - logLike - logPrior

    if np.log(rand()) < lap:
        x = px

    # preturb circle,  y
    px = x
    py = y + param_dict['yTau'] * randn()
    pz = z
    pr = r
    pxE = xE
    pyE = yE
    pzE = zE

    logLike = metrics.calc_log_like_sphere_mix(x, y, z, xE, yE, zE, r,  data[mles_param_dict['sInGIndex'], :],  mles_param_dict['rSig'])
    pLogLike = metrics.calc_log_like_sphere_mix(px, py, pz, pxE, pyE, pzE, pr,  data[mles_param_dict['sInGIndex'], :],  mles_param_dict['rSig'])
    logPrior = metrics.calc_log_prior_sphere(x, y, z, r, xE, yE, zE, prior)
    pLogPrior = metrics.calc_log_prior_sphere(px, py, pz, pxE, pyE, pzE, pr, prior)

    lap = pLogLike + pLogPrior - logLike - logPrior

    if np.log(rand()) < lap:
        y = py

    # preturb circle,  z
    px = x
    py = y
    pz = z + param_dict['zTau'] * randn()
    pr = r
    pxE = xE
    pyE = yE
    pzE = zE

    logLike = metrics.calc_log_like_sphere_mix(x, y, z, xE, yE, zE, r, data[mles_param_dict['sInGIndex'], :],  mles_param_dict['rSig'])
    pLogLike = metrics.calc_log_like_sphere_mix(px, py, pz, pxE, pyE, pzE, pr, data[mles_param_dict['sInGIndex'], :],  mles_param_dict['rSig'])
    logPrior = metrics.calc_log_prior_sphere(x, y, z, r, xE, yE, zE, prior)
    pLogPrior = metrics.calc_log_prior_sphere(px, py, pz, pxE, pyE, pzE, pr, prior)

    lap = pLogLike + pLogPrior - logLike - logPrior

    if np.log(rand()) < lap:
        z = pz

    # preturb axis,  x
    px = x
    py = y
    pz = z
    pr = r
    pxE = xE + param_dict['xTau'] * (1 - mles_param_dict['fix_xE'])*randn()
    while pxE <= 0:
        pxE = xE + param_dict['xETau'] * (1 - mles_param_dict['fix_xE'])*randn()

    pyE = yE
    pzE = zE

    logLike = metrics.calc_log_like_sphere_mix(x, y, z, xE, yE, zE, r, data[mles_param_dict['sInGIndex'], :],  mles_param_dict['rSig'])
    pLogLike = metrics.calc_log_like_sphere_mix(px, py, pz, pxE, pyE, pzE, pr, data[mles_param_dict['sInGIndex'], :],  mles_param_dict['rSig'])
    logPrior = metrics.calc_log_prior_sphere(x, y, z, r, xE, yE, zE, prior)
    pLogPrior = metrics.calc_log_prior_sphere(px, py, pz, pxE, pyE, pzE, pr, prior)

    lap = pLogLike + pLogPrior - logLike - logPrior

    if np.log(rand()) < lap:
        xE = pxE

    # preturb circle,  y
    px = x
    py = y
    pz = z
    pr = r
    pxE = xE
    pyE = yE + param_dict['yTau'] * (1 - mles_param_dict['fix_yE'])*randn()
    while (pyE <= 0):
        pyE = yE + param_dict['yTau'] * (1 - mles_param_dict['fix_yE'])*randn()

    pzE = zE

    logLike = metrics.calc_log_like_sphere_mix(x, y, z, xE, yE, zE, r, data[mles_param_dict['sInGIndex'], :],  mles_param_dict['rSig'])
    pLogLike = metrics.calc_log_like_sphere_mix(px, py, pz, pxE, pyE, pzE, pr, data[mles_param_dict['sInGIndex'], :],  mles_param_dict['rSig'])
    logPrior = metrics.calc_log_prior_sphere(x, y, z, r, xE, yE, zE, prior)
    pLogPrior = metrics.calc_log_prior_sphere(px, py, pz, pxE, pyE, pzE, pr, prior)

    lap = pLogLike + pLogPrior - logLike - logPrior

    if np.log(rand()) < lap:
        yE = pyE

    # preturb circle,  z
    px = x
    py = y
    pz = z
    pr = r
    pxE = xE
    pyE = yE
    pzE = zE + param_dict['zTau'] * (1 - mles_param_dict['fix_zE'])*randn()
    while (pzE <= 0):
        pzE = zE + param_dict['zTau'] * (1 - mles_param_dict['fix_zE'])*randn()

    logLike = metrics.calc_log_like_sphere_mix(x, y, z, xE, yE, zE, r, data[mles_param_dict['sInGIndex'], :],  mles_param_dict['rSig'])
    pLogLike = metrics.calc_log_like_sphere_mix(px, py, pz, pxE, pyE, pzE, pr, data[mles_param_dict['sInGIndex'], :],  mles_param_dict['rSig'])
    logPrior = metrics.calc_log_prior_sphere(x, y, z, r, xE, yE, zE, prior)
    pLogPrior = metrics.calc_log_prior_sphere(px, py, pz, pxE, pyE, pzE, pr, prior)

    lap = pLogLike + pLogPrior - logLike - logPrior

    if np.log(rand()) < lap:
        zE = pzE

    # preturb circle,  r
    px = x
    py = y
    pz = z
    pr = r + param_dict['rTau'] * (1 - mles_param_dict['fix_r'])*randn()
    while(pr <= 0):
        pr = r + param_dict['rTau'] * (1 - mles_param_dict['fix_r'])*randn()

    pxE = xE
    pyE = yE
    pzE = zE

    logLike = metrics.calc_log_like_sphere_mix(x, y, z, xE, yE, zE, r, data[mles_param_dict['sInGIndex'], :],  mles_param_dict['rSig'])
    pLogLike = metrics.calc_log_like_sphere_mix(px, py, pz, pxE, pyE, pzE, pr, data[mles_param_dict['sInGIndex'], :],  mles_param_dict['rSig'])
    logPrior = metrics.calc_log_prior_sphere(x, y, z, r, xE, yE, zE, prior)
    pLogPrior = metrics.calc_log_prior_sphere(px, py, pz, pxE, pyE, pzE, pr, prior)

    lap = pLogLike + pLogPrior - logLike - logPrior

    if np.log(rand()) < lap:
        r = pr

    # preturb full circle
    px = x + param_dict['xTau'] * randn()
    py = y + param_dict['yTau'] * randn()
    pz = z + param_dict['zTau'] * randn()
    pr = r + param_dict['rTau'] * (1 - mles_param_dict['fix_r'])*randn()
    while(pr <= 0):
        pr = r + param_dict['rTau'] * (1 - mles_param_dict['fix_r'])*randn()

    pxE = xE + param_dict['xTau'] * (1 - mles_param_dict['fix_xE'])*randn()
    while (pxE <= 0):
        pxE = xE + param_dict['xTau'] * (1 - mles_param_dict['fix_xE'])*randn()

    pyE = yE + param_dict['yTau'] * (1 - mles_param_dict['fix_yE'])*randn()
    while (pyE <= 0):
        pyE = yE + param_dict['yTau'] * (1 - mles_param_dict['fix_yE'])*randn()

    pzE = zE + param_dict['zTau'] * (1 - mles_param_dict['fix_zE'])*randn()
    while (pzE <= 0):
        pzE = zE + param_dict['zTau'] * (1 - mles_param_dict['fix_zE'])*randn()

    logLike = metrics.calc_log_like_sphere_mix(x, y, z, xE, yE, zE, r, data[mles_param_dict['sInGIndex'], :],  mles_param_dict['rSig'])
    pLogLike = metrics.calc_log_like_sphere_mix(px, py, pz, pxE, pyE, pzE, pr, data[mles_param_dict['sInGIndex'], :],  mles_param_dict['rSig'])
    logPrior = metrics.calc_log_prior_sphere(x, y, z, r, xE, yE, zE, prior)
    pLogPrior = metrics.calc_log_prior_sphere(px, py, pz, pxE, pyE, pzE, pr, prior)

    lap = pLogLike + pLogPrior - logLike - logPrior

    if np.log(rand()) < lap:
        x = px
        y = py
        z = pz
        r = pr
        xE = pxE
        yE = pyE
        zE = pzE

    return (x,  y,  z,  r,  xE,  yE,  zE)
