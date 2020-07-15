#!/usr/bin/env python

import json
import numpy as np
import os


class ParamSS:
    def __init__(self,  n_data,  meanR):
        self.xc = [0, 0]
        self.yc = [0, 0]
        self.zc = [0, 0]
        self.xE = [0, 0]
        self.yE = [0, 0]
        self.zE = [0, 0]
        self.r = [0, 0]
        self.rSig = [0, 0]
        self.sInG = np.zeros((n_data, 2))
        self.meanR = meanR

    def set_params(self,  param_dict):
        self.xc[0] = self.xc[0] + param_dict['xc']
        self.xc[1] = self.xc[1] + param_dict['xc']**2
        self.yc[0] = self.yc[0] + param_dict['yc']
        self.yc[1] = self.yc[1] + param_dict['yc']**2
        self.zc[0] = self.zc[0] + param_dict['zc']
        self.zc[1] = self.zc[1] + param_dict['zc']**2
        self.xE[0] = self.xE[0] + param_dict['xE']
        self.xE[1] = self.xE[1] + param_dict['xE']**2
        self.yE[0] = self.yE[0] + param_dict['yE']
        self.yE[1] = self.yE[1] + param_dict['yE']**2
        self.zE[0] = self.zE[0] + param_dict['zE']
        self.zE[1] = self.zE[1] + param_dict['zE']**2
        self.r[0] = self.r[0] + param_dict['r']
        self.r[1] = self.r[1] + param_dict['r']**2
        self.rSig[0] = self.rSig[0] + param_dict['rSig']
        self.rSig[1] = self.rSig[1] + param_dict['rSig']**2
        self.sInG[:, 0] = self.sInG[:, 0] + param_dict['sInG']
        self.sInG[:, 1] = self.sInG[:, 1] + param_dict['sInG']**2

    def summarize_params(self,  sample):
        self.xc[0] = self.xc[0]/sample
        self.xc[1] = np.sqrt(self.xc[1]/sample - self.xc[0]**2)
        self.yc[0] = self.yc[0]/sample
        self.yc[1] = np.sqrt(self.yc[1]/sample - self.yc[0]**2)
        self.zc[0] = self.zc[0]/sample
        self.zc[1] = np.sqrt(self.zc[1]/sample - self.zc[0]**2)
        self.xE[0] = self.xE[0]/sample
        self.xE[1] = np.sqrt(self.xE[1]/sample - self.xE[0]**2)
        self.yE[0] = self.yE[0]/sample
        self.yE[1] = np.sqrt(self.yE[1]/sample - self.yE[0]**2)
        self.zE[0] = self.zE[0]/sample
        self.zE[1] = np.sqrt(self.zE[1]/sample - self.zE[0]**2)
        self.r[0] = self.r[0]/sample
        self.r[1] = np.sqrt(self.r[1]/sample - self.r[0]**2)
        self.rSig[0] = self.rSig[0]/sample
        self.rSig[1] = np.sqrt(self.rSig[1]/sample - self.rSig[0]**2)
        self.sInG[:, 0] = self.sInG[:, 0]/sample
        self.sInG[:, 1] = np.sqrt(self.sInG[:, 1]/sample - self.sInG[:, 0]**2)

        self.sInGIndex = np.where(self.sInG[:, 0] > 0.5)[0]
        self.sOutGIndex = np.where(self.sInG[:, 0] <= 0.5)[0]

    def output_csv(self, output_dir,  cell_id, sep='\t'):
        base_file_name = '{}.csv'.format(cell_id)
        file_path = os.path.join(output_dir, base_file_name)
        with open(file_path,  'w') as fh:
            fh.write('Name{}X{}Y{}Z{}Scale{}rSig\n'.format(sep, sep, sep, sep, sep))
            fh.write('{}{}{:.6f}{}{:.6f}{}{:.6f}{}{:.6f}{}{:.6f}\n'.format(cell_id, sep, self.xc[0], sep, self.yc[0], sep, self.zc[0], sep, self.r[0], sep, self.rSig[0]))

    def output_json(self, output_dir, cell_id):
        base_file_name = '{}.json'.format(cell_id)
        file_path = os.path.join(output_dir, base_file_name)
        save_dict = self.__dict__
        save_dict['sInG'] = save_dict['sInG'].tolist()
        save_dict['sInGIndex'] = save_dict['sInGIndex'].tolist()
        save_dict['sOutGIndex'] = save_dict['sOutGIndex'].tolist()
        save_dict = {'paramSS': save_dict}
        with open(file_path, 'w') as fh:
            fh.write(json.dumps(save_dict))
