#!/usr/bin/env python

from data_utils import map_to_non_euclid_3d


def calc_log_like_sphere_mix(x, y, z, xE, yE, zE, r, data, rSig):
    # n = data.shape[0]
    # arcLengthSig = max(MixThetas)/2 # 2 Std deveations between points on arch
    rTheta = map_to_non_euclid_3d(data, x, y, z, xE, yE, zE)

    ll = - (0.5/(rSig*rSig))*sum((rTheta[:, 2]-r)**2)
    # Just look at radius from concentric circles ...
    # for i in range(n):
    #    tTheta = min(abs(rTheta(i, 2)-MixThetas))
    #    ll = ll - (0.5/(arcLengthSig*arcLengthSig))*(tTheta**2)
    return ll


def calc_log_prior_sphere(x, y, z, r, xE, yE, zE, prior):
    lp = -(0.5/(prior.xStd*prior.xStd))*((x-prior.xMean)**2)
    lp = lp - (0.5/(prior.yStd*prior.yStd))*((y-prior.yMean)**2)
    lp = lp - (0.5/(prior.zStd*prior.zStd))*((z-prior.zMean)**2)
    lp = lp - (0.5/(prior.rStd*prior.rStd))*((r-prior.rMean)**2)
    lp = lp - (0.5/(prior.xEStd*prior.xEStd))*((xE-prior.xEMean)**2)
    lp = lp - (0.5/(prior.yEStd*prior.yEStd))*((yE-prior.yEMean)**2)
    lp = lp - (0.5/(prior.zEStd*prior.zEStd))*((zE-prior.zEMean)**2)
    return lp
