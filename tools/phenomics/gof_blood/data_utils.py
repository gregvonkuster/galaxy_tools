#!/usr/bin/env python

import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt


def load_raw3d_data(filename):
    with open(filename, 'r') as file:
        data = np.loadtxt(file, delimiter=',', usecols=(0, 1, 2), skiprows=1)
    return data


def normalize_data(data, log_fh):
    dataNorm = np.copy(data)
    dataNorm[:, 0] = data[:, 0]-np.mean(data[:, 0])
    dataNorm[:, 1] = data[:, 1]-np.mean(data[:, 1])
    dataNorm[:, 2] = data[:, 2]-np.mean(data[:, 2])
    log_fh.write("data-mean\n{}\n".format(dataNorm[:5]))
    pca = PCA()
    dataNorm1 = pca.fit_transform(dataNorm)
    log_fh.write("data-pca\n{}\n".format(dataNorm1[:5]))

    # maxAll = np.amax(dataNorm1)*1.1
    # minAll = np.amin(dataNorm1)*1.1

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # X, Y, Z = dataNorm1[:, 0], dataNorm1[:, 1], dataNorm1[:, 2]
    ax.scatter(dataNorm1[:, 0], dataNorm1[:, 1], dataNorm1[:, 2])
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    plt.show(block=False)

    r = np.sqrt((dataNorm[:, 0]**2)+(dataNorm[:, 1]**2)+(dataNorm[:, 2]**2))
    meanR = np.mean(r)
    # stdR = np.std(r)
    dataNorm = dataNorm/meanR

    r = np.sqrt((dataNorm1[:, 0]**2)+(dataNorm1[:, 1]**2)+(dataNorm1[:, 2]**2))
    meanR = np.mean(r)
    # stdR = np.std(r)
    dataNorm1 = dataNorm1/meanR
    # maxAll = np.amax(dataNorm1)*1.1
    # minAll = np.amin(dataNorm1)*1.1
    plot_2D_graphs(dataNorm1)

    return dataNorm1, meanR


def map_to_non_euclid_3d(data, xc=0, yc=0, zc=0, xE=1, yE=1, zE=1):
    """
    Converts data form Catesian coordinate system to spherical coordinate systems
    data (numpy array): Set of points in 3D Cartesian Coordinates.
    Output (rTheta) - Data points in spherical coordinate systems
    """
    n = data.shape[0]
    rTheta = np.zeros(data.shape)
    rTheta[:, 2] = np.sqrt(((data[:, 0] - xc) / xE) ** 2 + ((data[:, 1] - yc) / yE) ** 2 + ((data[:, 2] - zc) / zE) ** 2)
    # Phi: acos(z/r)
    rTheta[:, 1] = np.arccos(np.divide(((data[:, 2] - zc) / zE), rTheta[:, 2]))
    # Theta: atan(y/x)
    rTheta[:, 0] = np.arctan(np.divide((abs(data[:, 1] - yc) / yE), ((data[:, 0] - xc) / xE)))

    for i in range(n):
        # Quadrant I doesn't change.

        if data[i, 0] < 0:
            if data[i, 1] > 0:
                # Quadrant II
                rTheta[i, 0] = np.pi + rTheta[i, 0]
            elif data[i, 1] < 0:
                # Quadrant III
                rTheta[i, 0] = np.pi - rTheta[i, 0]
        elif data[i, 0] > 0 and data[i, 1] < 0:
            # Quadrant IV
            rTheta[i, 0] = 2 * np.pi - rTheta[i, 0]

    return rTheta


def map_to_non_euclid_to_plot_3d(rTheta):
    tTheta = np.zeros(rTheta.shape)
    tTheta[:, 2] = rTheta[:, 2]
    tTheta[:, 1] = (rTheta[:, 1]*rTheta[:, 2])*np.sin(rTheta[:, 0])
    tTheta[:, 0] = (rTheta[:, 1]*rTheta[:, 2])*np.cos(rTheta[:, 0])
    return tTheta


def plot_2D_graphs(data):
    fig, axs = plt.subplots(3)
    axs[0].scatter(data[:, 0], data[:, 1], linestyle='--')
    axs[0].set_xlabel('X Axis')
    axs[0].set_ylabel('Y Axis')
    axs[0].grid()

    axs[1].scatter(data[:, 0], data[:, 2], color='y', linestyle='--')
    axs[1].set_xlabel('X Axis')
    axs[1].set_ylabel('Z Axis')
    axs[1].grid()

    axs[2].scatter(data[:, 1], data[:, 2], color='green', linestyle='--')
    axs[2].set_xlabel('Y Axis')
    axs[2].set_ylabel('Z Axis')
    axs[2].grid()
    plt.show(block=False)
