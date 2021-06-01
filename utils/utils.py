import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from sklearn.preprocessing import normalize

def blackbox_function(x, y=None, sim=False):
    if sim:
        if y is None:
            return -x ** 2 + 6
        else:
            return -(x+y) ** 2 + 6

    # Reading the magnitude of the N170 data
    filename = 'Output.txt'
    lines = open(filename).read().splitlines()

    try:
        latency = float(lines[-1])
    except ValueError:
        print('Failed to convert value to float')
        wait = input("PRESS ENTER TO CONTINUE.")
        lines = open(filename).read().splitlines()
        latency = float(lines[-1])
    except IndexError:
        print('The latent file is empty')
        wait = input("PRESS ENTER TO CONTINUE.")
        lines = open(filename).read().splitlines()
        latency = float(lines[-1])
    return latency


def obtain_confidence(sim=False):
    if sim:
        noise = np.random.normal(0, 0.60, size=1)[0]
        return noise

    # Reading the Confidence levels of the target value
    filename = 'Confidence.txt'
    lines = open(filename).read().splitlines()

    try:
        confidence = float(lines[-1])
    except ValueError:
        print('Failed to convert confidence value to float')
        wait = input("PRESS ENTER TO CONTINUE.")
        lines = open(filename).read().splitlines()
        confidence = float(lines[-1])
    except IndexError:
        print('The confidence file is empty')
        wait = input("PRESS ENTER TO CONTINUE.")
        lines = open(filename).read().splitlines()
        confidence = float(lines[-1])
    return confidence



def posterior(optimizer, x_obs, y_obs, grid):

    optimizer._gp.fit(x_obs, y_obs)

    mu, sigma = optimizer._gp.predict(grid, return_std=True)
    return mu, sigma


def plot_gp(optimizer, logpath, i, utility_function, bounds, x, y=None):

    fig = plt.figure(figsize=(16, 10))
    steps = len(optimizer.res)
    fig.suptitle(
        'Gaussian Process and Utility Function After {} Steps'.format(steps),
        fontdict={'size': 30}
    )

    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    axis = plt.subplot(gs[0])
    acq = plt.subplot(gs[1])

    # x_obs = np.array([[res["params"]["x"]] for res in optimizer.res])
    # y_obs = np.array([res["target"] for res in optimizer.res])

    x_obs = np.array([[res["params"]["x"]] for res in optimizer.res])
    y_obs = np.array([res["target"] for res in optimizer.res])

    y_obs, norm = normalize(y_obs.reshape(1, -1), return_norm=True)
    y_obs = y_obs.flatten()

    mu, sigma = posterior(optimizer, x_obs, y_obs, x)
    utility = utility_function.utility(x, optimizer._gp, y_obs.max())

    # Unnormalize data
    mu = mu*norm
    sigma = sigma*norm
    y_obs = y_obs*norm

    if y is not None:
        axis.plot(x, y, linewidth=3, label='Target')
    axis.plot(x_obs.flatten(), y_obs, 'D', markersize=8, label=u'Observations', color='r')
    axis.plot(x, mu, '--', color='k', label='Prediction')

    axis.fill(np.concatenate([x, x[::-1]]),
              np.concatenate([mu - 1.9600 * sigma, (mu + 1.9600 * sigma)[::-1]]),
              alpha=.6, fc='c', ec='None', label='95% confidence interval')
    # if(bounds == "large"):
    #     axis.set_xlim((-1, 1))
    # else:
    #     axis.set_xlim((0, 1))
    # axis.set_ylim((None, None))
    axis.set_ylabel('f(x)', fontdict={'size': 20})
    axis.set_xlabel('x', fontdict={'size': 20})

    # utility = utility_function.utility(x, optimizer._gp, 0)
    acq.plot(x, utility, label='Utility Function', color='purple')
    acq.plot(x[np.argmax(utility)], np.max(utility), '*', markersize=15,
             label=u'Next Best Guess', markerfacecolor='gold', markeredgecolor='k', markeredgewidth=1)

    # if (bounds == "large"):
    #     acq.set_xlim((-1, 1))
    # else:
    #     acq.set_xlim((0, 1))
    # acq.set_ylim((0, np.max(utility) + 0.5))
    acq.set_ylabel('Utility', fontdict={'size': 20})
    acq.set_xlabel('x', fontdict={'size': 20})

    axis.legend(loc=2, bbox_to_anchor=(1.01, 1), borderaxespad=0.)
    acq.legend(loc=2, bbox_to_anchor=(1.01, 1), borderaxespad=0.)

    fig.savefig(logpath+'/fig_{}'.format(i))

