import numpy as np

from data_utils import map_to_non_euclid_3d


def gen_selection_in_g_3d(data, mles_param_dict, prior):

	n = data.shape[0]
	sInG = np.zeros(n,)
	#lnUniformProb = -0.5*log(2*np.pi*mles_param_dict['rSig']*mles_param_dict['rSig']) - 0.5*(prior.stdLevelReference)**2
	lnUniformProb = -0.5*(prior.stdLevelRefference)**2

	rTheta = map_to_non_euclid_3d(data, mles_param_dict['xc'], mles_param_dict['yc'], mles_param_dict['zc'], mles_param_dict['xE'], mles_param_dict['yE'], mles_param_dict['zE'])
	#lnRProbInG = -0.5*log(2*pi*mles_param_dict.rSig*mles_param_dict.rSig) - (0.5/(rSig*rSig))*((rTheta[:, 0]-mles_param_dict['r'])**2)
	lnRProbInG = -(0.5/(mles_param_dict['rSig']*mles_param_dict['rSig']))*((rTheta[:, 2]-mles_param_dict['r'])**2)

	for i in range(n):
	    pIn = 1/(1+np.exp(max(lnUniformProb-lnRProbInG[i], -200)))
	    if np.random.rand() < pIn:
	        sInG[i] = 1

	sInGIndex = np.where(sInG == 1)[0]
	sOutGIndex = np.where(sInG == 0)[0]

	return sInG, sInGIndex, sOutGIndex
