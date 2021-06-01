import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


def plot_bo(mu, sigma, x, y, show = True, save = False, i = 0, logpath='logs/'):

    Zmu = np.reshape(mu, (50, 50))
    Zsigma = np.reshape(sigma, (50, 50))

    conf0 = np.array(mu - 1.9600 * sigma).reshape(50, 50)
    conf1 = np.array(mu + 1.9600 * sigma).reshape(50, 50)

    fig = plt.figure(1, figsize=(23, 8))
    X0p, X1p = np.meshgrid(x, y)

    font_dict_title = {'fontsize': 25}
    font_dict_label = {'fontsize': 18}
    font_dict_label3 = {'fontsize': 15}

    ax0 = fig.add_subplot(121)
    fig0 = ax0.pcolormesh(X0p, X1p, Zmu)
    ax0.set_title('Gaussian Process Predicted Mean', fontdict=font_dict_title)
    ax0.set_xlabel('Component 1', fontdict=font_dict_label)
    ax0.set_ylabel('Component 2', fontdict=font_dict_label)
    fig.colorbar(fig0)

    ax1 = fig.add_subplot(122)
    fig1 = ax1.pcolormesh(X0p, X1p, Zsigma)
    ax1.set_title('Gaussian Process Variance', fontdict=font_dict_title)
    ax1.set_xlabel('Component 1', fontdict=font_dict_label)
    ax1.set_ylabel('Component 2', fontdict=font_dict_label)
    fig.colorbar(fig1)

    fig_3d = plt.figure(2, figsize=(23, 8))

    ax2 = fig_3d.add_subplot(221, projection='3d')
    fig2 = ax2.plot_surface(X0p, X1p, Zmu, label='prediction', cmap=cm.coolwarm)
    ax2.set_title('Gaussian Process Mean', fontdict=font_dict_title)
    ax2.set_xlabel('Component 1', fontdict=font_dict_label3)
    ax2.set_ylabel('Component 2', fontdict=font_dict_label3)
    ax2.set_zlabel('P. Mean', fontdict=font_dict_label3)

    ax3 = fig_3d.add_subplot(222, projection='3d')
    fig3 = ax3.plot_surface(X0p, X1p, Zsigma, cmap=cm.coolwarm)
    ax3.set_title('Gaussian Process Variance', fontdict=font_dict_title)
    ax3.set_xlabel('Component 1', fontdict=font_dict_label3)
    ax3.set_ylabel('Component 2', fontdict=font_dict_label3)
    ax3.set_zlabel('Variance', fontdict=font_dict_label3)

    ax4 = fig_3d.add_subplot(223, projection='3d')
    fig4 = ax4.plot_surface(X0p, X1p, conf0, label='confidence', alpha=0.3)
    fig4 = ax4.plot_surface(X0p, X1p, conf1, label='confidence', alpha=0.3)
    ax4.set_title('95% Confidence Interval', fontdict=font_dict_title)
    ax4.set_xlabel('Component 1', fontdict=font_dict_label3)
    ax4.set_ylabel('Component 2', fontdict=font_dict_label3)
    ax4.set_zlabel('P.Mean', fontdict=font_dict_label3)

    if show:
        plt.show()
    if save:
        fig.savefig(logpath+'/2d_{}'.format(i))
        fig_3d.savefig(logpath+'/3d_{}'.format(i))

    fig.clf()
    fig_3d.clf()

    return

