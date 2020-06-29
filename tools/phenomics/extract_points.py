#!/usr/bin/env python

import csv
import os
import vtk
from vtk.util.numpy_support import vtk_to_numpy

INPUT_DIR = 'input_dir'
OUTPUT_DIR = 'output_dir'


def read_vtk(file_path):
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(file_path)
    reader.Update()
    polydata = reader.GetOutput()
    points = polydata.GetPoints()
    array = points.GetData()
    numpy_nodes = vtk_to_numpy(array)
    return numpy_nodes


def get_base_file_name(file_path):
    base_file_name = os.path.basename(file_path)
    if base_file_name.find('.') > 0:
        # Eliminate the extension.
        return os.path.splitext(base_file_name)[0]
    return base_file_name


for file_name in sorted(os.listdir(INPUT_DIR)):
    file_path = os.path.abspath(os.path.join(INPUT_DIR, file_name))
    base_file_name = get_base_file_name(file_path)
    coordinate = read_vtk(file_path)
    if len(coordinate) == 0:
        pass
    else:
        out_file_path = os.path.abspath(os.path.join(OUTPUT_DIR, '%s.csv' % base_file_name))
        with open(out_file_path, 'w', newline='', encoding='utf-8-sig') as fh:
            csv_writer = csv.writer(fh)
            head = ['x' 'y' 'z']
            for l in head:
                csv_writer.writerow(l)
            for l in coordinate:
                csv_writer.writerow(l)
