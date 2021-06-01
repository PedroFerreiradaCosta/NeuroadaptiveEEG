# This script uses Bayesian Optimization to navigate in the face Space.
# It allows to add a figure (in folder 'new_face') in order to add them
# in real time to the space.

import numpy as np
import os
import sys
import glob
import time
import argparse
import joblib

from PIL import Image
from sklearn.gaussian_process.kernels import Matern, WhiteKernel
from utils.bayesian_optimization import Bayesian_Optimization, UtilityFunction

from utils.utils import plot_gp
if not os.path.exists('stimulus/'):
    os.makedirs('stimulus/')

parser = argparse.ArgumentParser(description="Inputs for the formatting of the Bayesian Optimizer.")
# parser.add_argument('-i', dest='init_points', type=int, nargs=1, default=[5],
#                     help="pass the number of random points sampled from the space.")
parser.add_argument('-n', dest='n_iter', type=int, nargs=1, default=[15],
                    help='pass the number of iterations for the Optimizer to run.')
parser.add_argument('-d', dest='face_display', type=bool, nargs='?', default=True,
                    help='Display faces sampled by the Optimizer.')
# parser.add_argument('-f', dest='save_file', type=bool, nargs='?', default=True,
#                     help='Name of Log file to save iteration data.')
parser.add_argument('-v', dest='vis', type=bool, nargs='?', default=True,
                    help='Show 2D and 3D visualization of the algorithm\'s progress.')
parser.add_argument('-s', dest='sim', type=bool, nargs='?',
                    help='Simulate target function.')

####################################################################################
# Example
# python3 1dBO_optm_signal.py -n 5 -d False -v True -s True
# python3 1dBO_optm_signal.py -s True
# python3 1dBO_optm_signal.py
#####################################################################################

def optm_signal(n_iter=[15], face_display=True, vis=False,
 sim=False):
    """
    This function runs the Bayesian Optimisation while it registers target values from file n170output.txt

    :type init_points: int
    :type n_iter:  int
    :type face_display: bool
    :type save_file: bool
    :type sim: bool
    :param init_points: number of points sampled randomly before the algorithm is run
    :param n_iter:  number of points samples by the algorithm
    :param face_display: if the image correspondent to the point sampled is displayed or not
    :param save_file: if a file with the iteration results is saved or not
    :param sim: if the target function is simulated or not.
    :return:
    """

    try: 
        n_iter = n_iter[0]
        np.random.seed(99)
        utility_function = 'ei'

        filename = 'Subj_ID.txt'
        lines = open(filename).read().splitlines()
        try:
            save_file = str(lines[-1])
        except IndexError:
            print("Failed to read participant name. Will save data in participant_x")
            save_file = 'participant_x'

        log_dir = 'data/'+ save_file + '/'
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
        log_name = '1dbo_' + str(int(round(time.time())))
        LOG_PATH = log_dir+log_name+'/'


        if face_display:
            # Change to Path?
            images_path = sorted(list(glob.glob('1dBO/images/*.png')))  # PATH TO 9 FACES

        space = np.array([-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2 ,0.4, 0.6, 0.8, 1.0, 1.2])

        if sim:
            def blackbox_function(x):
                return -x ** 2 + 1 + np.random.normal(0,0.1,1)[0]
        else:
            def blackbox_function(x):

                # Reading the magnitude of the N170 data
                filename = './Output.txt'
                lines = open(filename).read().splitlines()

                try:
                    latency = float(lines[-1])
                except ValueError:
                    print('Failed to convert value to float')
                    wait = input("PRESS ENTER TO CONTINUE.")
                    lines = open(filename).read().splitlines()
                    latency = float(lines[-1])
                except IndexError:
                    wait = input("PRESS ENTER TO CONTINUE.")
                    lines = open(filename).read().splitlines()
                    latency = float(lines[-1])

                return latency

        pbounds = {'x': (-1, 1.2)}
        optimizer = Bayesian_Optimization(
            pbounds=pbounds,
            random_state=42,
        )
        kernel = Matern(length_scale=0.2, length_scale_bounds=(1e-1, 1e1), nu=1.5) + WhiteKernel(noise_level=0.5)
        optimizer._gp.set_params(kernel=kernel, normalize_y=True)

        # 2 Utility functions w different acquisition functions w parameters to optimize exploration
        if utility_function == 'ei':
            utility = UtilityFunction(function="ei", hyperparam=1e-1) #0 to 0.1
        elif utility_function == 'ucb':
            # check hps
            utility = UtilityFunction(function="ucb", hyperparam=1) # 1 to 10 - bigger more exploratory

        prev_latency = 0
        warm_start = np.array([1.2, -1, 0.2, -0.2])
        init_points = warm_start.shape[0]

        if not os.path.exists(LOG_PATH):
            os.makedirs(LOG_PATH)

        with open(LOG_PATH + 'log.txt', "w+") as f:
            f.write('|   iter    |  target   |     x     |     y     |\n')
            f.write('-------------------------------------------------\n')

            for it in range(init_points + n_iter):
                if it < init_points:
                    next_point_to_probe = {'x': warm_start[it]}

                    print("Next point to probe is:", next_point_to_probe)

                    x = list(next_point_to_probe.values())
                else:
                    next_point_to_probe = optimizer.suggest(utility,
                                                            n_warmup=1000,
                                                            n_iter=100)

                    print("Next point to probe is:", next_point_to_probe)

                    x = list(next_point_to_probe.values())

                for i in range(0, 10):
                    #step = i * 0.01
                    step = i * 0.1
                    close_points_idx = np.where((space[:] > x[0] - step) & (space[:] < x[0] + step))
                    if close_points_idx[0].size > 0:  # Check if not empty
                        break

                # If more than one point identified, use euclidean distance to choose best
                a = np.array(x)
                dist = []
                idx_dist = 0

                for i in range(0, close_points_idx[0].size):
                    b = space[close_points_idx[0][i]]
                    tmp = np.linalg.norm(a - b)
                    dist.append(tmp)
                if len(dist) > 1:
                    idx_dist = np.where(dist == np.min(dist))
                print(f"Closest point: {space[close_points_idx[0][i]]}")

                idx_x = np.where(space[:] == space[close_points_idx[0][i]])[0]
                print(f"Number of sampled image: {idx_x}")
                # print('Real image position: {}'.format(b))
                print('Sampled position: {}'.format(a))
                # print('Distance: {}'.format(np.min(dist)))

                if face_display:
                    # identify image associated with closest point and save as new stimulus
                    try:
                        image_path = images_path[idx_x]
                    # when there is only 1 value, idx is returned as an array
                    except:
                        image_path = images_path[idx_x[0]]
                    img = Image.open(image_path)
                    img.save('stimulus/ganstim01.png')
                    img.show()


                target = blackbox_function(**next_point_to_probe)
                # 999 is the break signal from matlab script
                if target == 999:
                    break

                while target == prev_latency or target == 0:
                    time.sleep(1)
                    print('Waiting for update')
                    target = blackbox_function(**next_point_to_probe)
                prev_latency = target
                print("Found the target value to be: ", target)
                print("________________________________________")

                if sim:
                    # register_sample = {'x': b[0]}
                    register_sample = next_point_to_probe
                else:
                    register_sample = {'x': b} # The value registered is the existing point closest to the chosen
                try:
                    optimizer.register(
                    params=register_sample,
                    target=target,
                )
                except KeyError:
                    register_sample['x'].update(register_sample['x']+np.random.normal(0, 0.01, size=1)[0])
                    print("Data point is not unique")
                    optimizer.register(
                        params=register_sample,
                        target=target,
                    )

                if vis:
                    x = np.linspace(-1.1, 1.3, 50).reshape(-1, 1)
                    plot_gp(optimizer, LOG_PATH, it, utility, pbounds, x)


                if save_file:
                    f.write('|     {}     |   {:.3f}   |   {:.3f}   |\n'.format(
                        it, target, register_sample['x']))
                    f.write('----------------------------------------\n')

                if optimizer._params.shape[0] > 2:
                    if optimizer._params[-1] == optimizer._params[-2] and optimizer._params[-1] == optimizer._params[-3]:
                        break

            print("Maximum obtained by the algorithm:")
            print(optimizer.max)

            if save_file:
                f.write('Maximum obtained by the algorithm: x = {:.3f}'.format(optimizer.max['params']['x']))
        joblib.dump(optimizer, LOG_PATH+'optimizer.joblib')

    except KeyboardInterrupt as e:
        print(e)
        print("Maximum obtained by the algorithm:")
        print(optimizer.max)

        # f.write('Maximum obtained by the algorithm: x = {:.3f}'.format(optimizer.max['params']['x']))
        joblib.dump(optimizer, LOG_PATH+'optimizer.joblib')
        return
    return

if __name__ == '__main__':
    args = parser.parse_args(sys.argv[1:])
    optm_signal(n_iter=args.n_iter, face_display=args.face_display,
                vis=args.vis, sim=args.sim)
