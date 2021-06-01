import joblib
from sklearn.gaussian_process.kernels import Matern, WhiteKernel
from utils.bayesian_optimization import Bayesian_Optimization, UtilityFunction
from utils.utils import plot_gp, posterior
from sklearn.preprocessing import normalize
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = "Arial"


file_names = ['results/P_01/optimizer.joblib',
'results/P_02/optimizer.joblib',
'results/P_03/optimizer_no_habituation.joblib',
'results/P_04/optimizer.joblib',
'results/P_05/optimizer.joblib',
]



optimizer_ls = []
for file_name in file_names:
    optimizer_ls.append(joblib.load(file_name))

optimizer_habituation = joblib.load('results/P_03/optimizer.joblib')

utility = UtilityFunction(function="ei", hyperparam=1e-1) #0 to 0.1

x = np.linspace(-1.1, 1.3, 50).reshape(-1, 1)
pbounds = {'x': (-1, 1.2)}
kernel = Matern(length_scale=0.2, length_scale_bounds=(1e-1, 1e1), nu=1.5) + WhiteKernel(noise_level=0.5)
# optimizer._gp.set_params(kernel=kernel, normalize_y=True)
LOG_PATH = 'results/'

pos = [0]
neg = [1,2,3,4]
axis= plt.subplot(111)
# Habituation on P_03
x_obs = np.array([[res["params"]["x"]] for res in optimizer_habituation.res])
y_obs = np.array([res["target"] for res in optimizer_habituation.res])
axis.scatter(x_obs, -y_obs, c='black')
axis.set_title('Individual statistical models for participant P_03',
                 pad=50,
                fontdict={'size': 15, 'fontname':'Arial'},
                weight="bold")

color = ['#001fcc', '#4d67ff','#8093ff']

for i in range(3):

    x_obs = np.array([[res["params"]["x"]] for res in optimizer_habituation.res[:7+i]])
    y_obs = np.array([res["target"] for res in optimizer_habituation.res[:7+i]])
    axis.scatter(x_obs[-1], -y_obs[-1], c=color[i])
    # Normalize data
    y_obs, norm = normalize(y_obs.reshape(1, -1), return_norm=True)
    y_obs = y_obs.flatten()

    mu, sigma = posterior(optimizer_habituation, x_obs, y_obs, x)
    # Rescale data
    mu = mu*norm
    mu = -mu
    sigma = -sigma
    sigma = sigma*norm
    y_obs = y_obs*norm
    axis.plot(x, mu, c=color[i])

plt.legend(['After 7 iterations', '8th iteration', '9th iteration'])
right_side = axis.spines["right"]
right_side.set_visible(False)
right_side = axis.spines["top"]
right_side.set_visible(False)
plt.xticks([-1,0.1,1.2],['stranger','50% mother', 'mother'])
plt.show()
plt.savefig('results/image_paper2.png', transparent=True)


# Code to generate Figure 4

mu_ls = []
std_ls = []
y_obs_ls = []
x_obs_ls = []
for it in range(len(optimizer_ls)):

    x_obs = np.array([[res["params"]["x"]] for res in optimizer_ls[it].res])
    y_obs = np.array([res["target"] for res in optimizer_ls[it].res])

    # Normalize data
    y_obs, norm = normalize(y_obs.reshape(1, -1), return_norm=True)
    y_obs = y_obs.flatten()

    mu, sigma = posterior(optimizer_ls[it], x_obs, y_obs, x)

    # Rescale data
    mu = mu*norm
    if it in neg:
        mu = -mu
        sigma = -sigma
    sigma = sigma*norm
    y_obs = y_obs*norm

    mu_ls.append(mu)
    std_ls.append(sigma)
    y_obs_ls.append(y_obs)
    x_obs_ls.append(x_obs)


# plotting all data distinguishing by participant and polarity
participant = ['--', '-.',  '-',  '.',  ':']
# legend = ['BP_15', 'BP_16', 'BP_18', 'BP_18', 'BP_20', 'BP_20', 'BP_21', 'BP_21']
legend = ['P_01', 'P_02', 'P_03', 'P_04', 'P_05']
polarity = ['red', 'blue', 'blue', 'blue', 'blue']

axis = plt.subplot(111)

new_legend = []
for it in range(len(file_names)):
    axis.plot(x, mu_ls[it], participant[it], c=polarity[it])
    new_legend.append(legend[it])

    plt.fill(np.concatenate([x, x[::-1]]),
          np.concatenate([mu_ls[it] -  std_ls[it]/y_obs_ls[it].shape[0], (mu_ls[it] + std_ls[it]/y_obs_ls[it].shape[0])[::-1]]),
          alpha=.1, c=polarity[it], fc='k', ec='None', label='95% confidence interval')
plt.legend(new_legend)
plt.xticks([-1,0.1,1.2],['stranger','50% mother', 'mother'])

plt.legend(new_legend)
plt.xticks([-1,0.1,1.2],['stranger','50% mother', 'mother'])
axis.set_title('Individual statistical models to the participants target response',
                 pad=30,
                fontdict={'size': 15, 'fontname':'Arial'},
                weight="bold")
plt.ylabel('Nc Amplitude ($\\mu$V)', fontdict={'fontname':'Arial'})
right_side = axis.spines["right"]
right_side.set_visible(False)
right_side = axis.spines["top"]
right_side.set_visible(False)

plt.show()




# Creating Figure 4
from matplotlib import gridspec
y=None
hps = [1e-1, 5e-1, 5e0]
x_obs = np.array([[res["params"]["x"]] for res in optimizer_ls[2].res[:4]])
y_obs = np.array([res["target"] for res in optimizer_ls[2].res[:4]])

fig = plt.figure(figsize=(16, 10))

gs = gridspec.GridSpec(3, 3)#, height_ratios=[3,1])
axis = plt.subplot(gs[:,0])
acq1 = plt.subplot(gs[0,1])
acq2 = plt.subplot(gs[1,1])
acq3 = plt.subplot(gs[2,1])
hab = plt.subplot(gs[:,2])
acq = [acq1, acq2, acq3]

# Participant P_03
x_obs = np.array([[res["params"]["x"]] for res in optimizer_ls[2].res[:4]])
y_obs = np.array([res["target"] for res in optimizer_ls[2].res[:4]])

y_obs, norm = normalize(y_obs.reshape(1, -1), return_norm=True)
y_obs = y_obs.flatten()

mu, sigma = posterior(optimizer_ls[2], x_obs, y_obs, x)

# Unnormalize data
mu = mu*norm
sigma = sigma*norm
y_obs = y_obs*norm
star_color = ['silver', 'green', 'gold']
if y is not None:
    axis.plot(x, y, linewidth=3, label='Target')
axis.plot(x_obs.flatten(), -y_obs, '.', markersize=15, label=u'Observations', color='k')
axis.plot(x, -mu, '--', color='blue', label='Prediction')

axis.fill(np.concatenate([x, x[::-1]]),
          np.concatenate([-mu - 1.9600 * sigma, (-mu + 1.9600 * sigma)[::-1]]),
          alpha=.1, c='blue', fc='c', ec='None', label='95% confidence interval')

axis.set_ylabel('Nc Amplitude ($\\mu$V)', fontdict={'size': 12})
axis.set_xticks([])
axis.set_xticks([-1,0.1,1.2])
axis.set_xticklabels(['stranger','50% mother', 'mother'])
# axis.set_xlabel('x', fontdict={'size': 20})
axis.set_title('Surrogate Model',pad=30,
                fontdict={'size': 15, 'fontname':'Arial'},
                weight="bold")
acq[0].set_title('Acquisition function',pad=30,
                fontdict={'size': 15, 'fontname':'Arial'},
                weight="bold")
for i, hp in enumerate(hps):

    utility_function = UtilityFunction(function="ei", hyperparam=hp) #0 to 0.1

    # x_obs = np.array([[res["params"]["x"]] for res in optimizer.res])
    # y_obs = np.array([res["target"] for res in optimizer.res])

    utility = utility_function.utility(x, optimizer_ls[2]._gp, y_obs.max()/norm)

    acq[i].plot(x, utility, label='Utility Function', color='blue')
    acq[i].plot(x[np.argmax(utility)], np.max(utility), '*', markersize=12,
             label=u'Next Point to sample', markerfacecolor=star_color[i], markeredgecolor='k', markeredgewidth=1)
    axis.plot(x[np.argmax(utility)], -mu[np.argmax(utility)], '*', markersize=12,
             label=f'Next sample for $\\xi$ = {hp}', markerfacecolor=star_color[i], markeredgecolor='k', markeredgewidth=1)
    acq[i].set_ylabel(f'$\\xi$ = {hp}', fontdict={'size': 12})
    # acq[i].set_xlabel('x', fontdict={'size': 20})
    acq[i].set_xticks([-1,0.1,1.2])
    # right_side = acq[i].spines["right"]
    # right_side.set_visible(False)
    # right_side = acq[i].spines["top"]
    # right_side.set_visible(False)
    # right_side = acq[i].spines["left"]
    # right_side.set_visible(False)
    acq[i].set_xticklabels(['stranger','50% mother', 'mother'])
    acq[i].set_yticks([])
    acq[i].legend()
    right_side = acq[i].spines["right"]
    right_side.set_visible(False)
    right_side = acq[i].spines["top"]
    right_side.set_visible(False)
    right_side = acq[i].spines["left"]
    right_side.set_visible(False)


axis.legend(loc="lower left")#loc=2, bbox_to_anchor=(1.01, 1), borderaxespad=0.)
right_side = axis.spines["right"]
right_side.set_visible(False)
right_side = axis.spines["top"]
right_side.set_visible(False)

# Habituation on participant 18
x_obs = np.array([[res["params"]["x"]] for res in optimizer_habituation.res[:7]])
y_obs = np.array([res["target"] for res in optimizer_habituation.res[:7]])
hab.scatter(x_obs, -y_obs, c='black')
hab.set_title('Individual statistical models \n for participant P_03',
                 pad=50,
                fontdict={'size': 15, 'fontname':'Arial'},
                weight="bold")
color = ['#001fcc', '#4d67ff','#8093ff']
it = 2 # BP_18
for i in range(3):
    print(file_names[it])
    x_obs = np.array([[res["params"]["x"]] for res in optimizer_habituation.res[:7+i]])
    y_obs = np.array([res["target"] for res in optimizer_habituation.res[:7+i]])
    hab.scatter(x_obs[-1], -y_obs[-1], c=color[i])
    # Normalize data
    y_obs, norm = normalize(y_obs.reshape(1, -1), return_norm=True)
    y_obs = y_obs.flatten()

    mu, sigma = posterior(optimizer_habituation, x_obs, y_obs, x)
    # Rescale data
    mu = mu*norm
    mu = -mu
    sigma = -sigma
    sigma = sigma*norm
    y_obs = y_obs*norm
    hab.plot(x, mu, c=color[i])

hab.set_ylabel('Nc Amplitude ($\\mu$V)', fontdict={'fontname':'Arial'})
hab.legend(['After 7 iterations', '8th iteration', '9th iteration'])
right_side = hab.spines["right"]
right_side.set_visible(False)
right_side = hab.spines["top"]
right_side.set_visible(False)
hab.set_xticks([-1,0.1,1.2],['stranger','50% mother', 'mother'])
hab.set_ytitle

plt.autoscale(True)
plt.show()

plt.savefig('results/image_paper4.png', transparent=True)
