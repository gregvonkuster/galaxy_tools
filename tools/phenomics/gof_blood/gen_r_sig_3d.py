import numpy as np
from numpy.random import gamma
import time

from data_utils import map_to_non_euclid_3d

def gen_r_sig_3d(data, mles_param_dict, prior):
	rTheta = map_to_non_euclid_3d(data[mles_param_dict['sInGIndex'], :], mles_param_dict['xc'], mles_param_dict['yc'], mles_param_dict['zc'], mles_param_dict['xE'], mles_param_dict['yE'], mles_param_dict['zE'])

	n = mles_param_dict['sInGIndex'].shape[0]
	sh = 0.5*n + prior.rSigShape
	sc = 0.5*np.sum((rTheta[:, 2]-mles_param_dict['r'])**2) + prior.rSigScale

	rSig = np.sqrt(1/gamma(sh, 1/sc))
	return rSig

dur = 4000

# print(time.strftime('%H:%M:%S', time.gmtime(dur)))

