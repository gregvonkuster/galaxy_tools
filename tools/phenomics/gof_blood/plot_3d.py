import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from data_utils import map_to_non_euclid_3d,  map_to_non_euclid_to_plot_3d


def plot_3d(prjCtrl, data, param):
    #  Generate and plot longitudinal circles (pole to pole)
    longitudeC = np.zeros((prjCtrl.nSpineDots, prjCtrl.nSpineCircles, 3))
    latitudeC = np.zeros((prjCtrl.nSpineDots, prjCtrl.nSpineCircles, 3))
    for i in range(prjCtrl.nSpineCircles):
        #  x longitude
        longitudeC[:, i, 0] = param['r']*param['xE']*(np.sin(prjCtrl.longitudeC)*np.cos(prjCtrl.MixPhiMean[i])) + param['xc']
        #  y longitude
        longitudeC[:, i, 1] = param['r']*param['yE']*(np.sin(prjCtrl.longitudeC)*np.sin(prjCtrl.MixPhiMean[i])) + param['yc']
        #  z longitude
        longitudeC[:, i, 2] = param['r']*param['zE']*np.cos(prjCtrl.longitudeC) + param['zc']
        #  x longitude
        latitudeC[:, i, 0] = param['r']*param['xE']*(np.sin(prjCtrl.MixPhiMean[i])*np.cos(prjCtrl.longitudeC)) + param['xc']
        #  y longitude
        latitudeC[:, i, 1] = param['r']*param['yE']*(np.sin(prjCtrl.MixPhiMean[i])*np.sin(prjCtrl.longitudeC)) + param['yc']
        #  z longitude
        latitudeC[:, i, 2] = param['r']*param['zE']*np.cos(prjCtrl.MixPhiMean[i]) + param['zc']

    fig = plt.figure()
    ax = fig.add_subplot(111,  projection='3d')
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    for i in range(prjCtrl.nSpineCircles):
        ax.scatter(latitudeC[:, i, 0], latitudeC[:, i, 1], latitudeC[:, i, 2])
        # plt.show(block=False)
        # plot3(latitudeC[:, i, 0], latitudeC[:, i, 1], latitudeC[:, i, 2], 'r.')
        # if i==1:
        #   # plt.hold()
        # plot3(longitudeC[:, i, 0], longitudeC[:, i, 1], longitudeC[:, i, 2], 'r.')
        ax.scatter(longitudeC[:, i, 0], longitudeC[:, i, 1], longitudeC[:, i, 2])
        # plt.show(block=False)

    ax.scatter(data[param['sInGIndex'], 0], data[param['sInGIndex'], 1], data[param['sInGIndex'], 2])
    # plt.show(block=False)
    ax.scatter(data[param['sOutGIndex'], 0], data[param['sOutGIndex'], 1], data[param['sOutGIndex'], 2])
    # maxAll = max(max(data))*1.1
    # minAll = min(min(data))*1.1
    # axis([minAll maxAll minAll maxAll minAll maxAll])
    plt.show(block=False)
    # plt.hold()

    fig = plt.figure()
    ax = fig.add_subplot(111,  projection='3d')
    for i in range(prjCtrl.nSpineCircles):
        tData = np.array([latitudeC[:, i, 0],  latitudeC[:, i, 1],  latitudeC[:, i, 2]])
        tData = map_to_non_euclid_3d(tData, param['xc'], param['yc'], param['zc'], param['xE'], param['yE'], param['zE'])
        tData = map_to_non_euclid_to_plot_3d(tData)
        ax.scatter(tData[:, 0], tData[:, 1], tData[:, 2])
        # if i==1:
        #   # plt.hold()
        tData = np.array([longitudeC[:, i, 0],  longitudeC[:, i, 1],  longitudeC[:, i, 2]])
        tData = map_to_non_euclid_3d(tData, param['xc'], param['yc'], param['zc'], param['xE'], param['yE'], param['zE'])
        tData = map_to_non_euclid_to_plot_3d(tData)
        ax.scatter(tData[:, 0], tData[:, 1], tData[:, 2])

    nonEuclid = map_to_non_euclid_3d(data, param['xc'], param['yc'], param['zc'], param['xE'], param['yE'], param['zE'])
    nonEuclidToPlot = map_to_non_euclid_to_plot_3d(nonEuclid)
    ax.scatter(nonEuclidToPlot[param['sInGIndex'], 0], nonEuclidToPlot[param['sInGIndex'], 1], nonEuclidToPlot[param['sInGIndex'], 2])
    ax.scatter(nonEuclidToPlot[param['sOutGIndex'], 0], nonEuclidToPlot[param['sOutGIndex'], 1], nonEuclidToPlot[param['sOutGIndex'], 2])
    # maxAll = max(max(nonEuclidToPlot))*1.1
    # minAll = min(min(nonEuclidToPlot))*1.1
    # maxR = max(nonEuclidToPlot[:, 2])*1.1
    # minR = min(nonEuclidToPlot[:, 2])*0.9
    # axis([minAll maxAll minAll maxAll minR maxR])
    plt.show(block=False)
    # plt.hold()

    # maxAll = max(max(data))*1.1
    # minAll = min(min(data))*1.1
    fig,  axs = plt.subplots(3)
    axs[0].scatter(data[param['sInGIndex'], 0], data[param['sInGIndex'], 1], marker='.', c='b')
    # plt.hold()
    axs[0].scatter(data[param['sOutGIndex'], 1], data[param['sOutGIndex'], 2], marker='x',  c='r')
    for i in range(prjCtrl.nSpineCircles):
        axs[0].scatter(latitudeC[:, i, 1],  latitudeC[:, i, 2],  marker='.', c='r')
        axs[0].scatter(longitudeC[:, i, 1],  longitudeC[:, i, 2],  marker='.', c='r')

    # plt.hold()
    axs[0].set_xlabel('X Axis')
    axs[0].set_ylabel('Y Axis')
    # axis([minAll maxAll minAll maxAll])
    # subplot(3, 1, 2)
    axs[1].scatter(data[param['sInGIndex'], 0], data[param['sInGIndex'], 2], marker='.', c='b')
    # plt.hold()
    axs[1].scatter(data[param['sOutGIndex'], 0], data[param['sOutGIndex'], 2], marker='x', c='r')
    for i in range(prjCtrl.nSpineCircles):
        axs[1].scatter(latitudeC[:, i, 0],  latitudeC[:, i, 2],  marker='.', c='r')
        axs[1].scatter(longitudeC[:, i, 0],  longitudeC[:, i, 2],  marker='.', c='r')

    # plt.hold()
    axs[1].set_xlabel('X Axis')
    axs[1].set_ylabel('Z Axis')
    # axis([minAll maxAll minAll maxAll])

    # subplot(3, 1, 3)
    axs[2].scatter(data[param['sInGIndex'], 1], data[param['sInGIndex'], 2], marker='.', c='b')
    # plt.hold()
    axs[2].scatter(data[param['sOutGIndex'], 1], data[param['sOutGIndex'], 2], marker='x', c='r')
    for i in range(prjCtrl.nSpineCircles):
        axs[2].scatter(latitudeC[:, i, 1],  latitudeC[:, i, 2],  marker='.', c='r')
        axs[2].scatter(longitudeC[:, i, 1],  longitudeC[:, i, 2],  marker='.', c='r')

    # plt.hold()
    axs[2].set_xlabel('Y Axis')
    axs[2].set_ylabel('Z Axis')
    plt.show(block=False)

    # axis([minAll maxAll minAll maxAll])
