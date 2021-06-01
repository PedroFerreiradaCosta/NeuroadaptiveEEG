from sklearn.gaussian_process.kernels import Matern
from sklearn.gaussian_process import GaussianProcessRegressor
import warnings
import numpy as np
from scipy.stats import norm
from scipy.optimize import minimize
from sklearn.preprocessing import normalize

class Bayesian_Optimization():
    def __init__(self, pbounds, random_state=42):
        self._random_state = np.random.RandomState(random_state)

        # Get the name of the parameters
        self._keys = sorted(pbounds)
        # Create an array with parameters bounds
        self._bounds = np.array(
            [item[1] for item in sorted(pbounds.items(), key=lambda x: x[0])],
            dtype=np.float
        )
        self.dim = len(pbounds)
        # preallocated memory for X and Y points
        self._params = np.empty(shape=(0, self.dim))
        self._target = np.empty(shape=(0))

        # Default gp
        self._gp = GaussianProcessRegressor(
            kernel=Matern(nu=2.5),
            alpha=1e-6,
            normalize_y=True,
            n_restarts_optimizer=25,
            random_state=self._random_state,
        )

    @property
    def max(self):
        try:
            res = {
                'target': self._target.max(),
                'params': dict(
                    zip(self._keys, self._params[self._target.argmax()])
                )
            }
        except ValueError:
            res = {}
        return res

    @property
    def res(self):
        """Get all target values found and corresponding parametes."""
        params = [dict(zip(self._keys, p)) for p in self._params]

        return [
            {"target": target, "params": param}
            for target, param in zip(self._target, params)
        ]

    def params_to_array(self, params):
        try:
            assert set(params) == set(self._keys)
        except AssertionError:
            raise ValueError(
                "Parameters' keys ({}) do ".format(sorted(params)) +
                "not match the expected set of keys ({}).".format(self.keys)
            )
        return np.asarray([params[key] for key in self._keys])

    def array_to_params(self, x):
        try:
            assert len(x) == len(self._keys)
        except AssertionError:
            raise ValueError(
                "Size of array ({}) is different than the ".format(len(x)) +
                "expected number of parameters ({}).".format(len(self.keys))
            )
        return dict(zip(self._keys, x))

    def register(self, params, target):
        """Expect observation with known target"""
        try:
            x = np.asarray(params, dtype=float)
        except TypeError:
            x = self.params_to_array(params)

        # if x in self._params:
        #     raise KeyError('Data point {} is not unique'.format(x))

        self._params = np.concatenate([self._params, x.reshape(1, -1)])
        self._target = np.concatenate([self._target, [target]])

    def suggest(self, utility_class, n_warmup=500, n_iter=10):
        """
        Find point to sample from utility function max.
        """

        # Sklearn's GP throws a large number of warnings at times, but
        # we don't really need to see them here.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # Normalize target values
            #
            _target_norm = normalize(self._target.reshape(1, -1))
            _target_norm = _target_norm.flatten()
            self._gp.fit(self._params, _target_norm)


        suggestion = self.acq_max(
            utility_func=utility_class.utility,
            gp=self._gp,
            y_max=_target_norm.max(),
            bounds=self._bounds,
            random_state=self._random_state,
            n_warmup=n_warmup,
            n_iter=n_iter
        )

        return self.array_to_params(suggestion)

    def acq_max(self, utility_func, gp, y_max, bounds,
                random_state, n_warmup, n_iter):
        # Warmup
        x_rand = random_state.uniform(bounds[:,0], bounds[:,1],
                                      size=(n_warmup, bounds.shape[0])
                                      )
        y_rand = utility_func(x_rand, gp=gp, y_max=y_max)
        x_max = x_rand[y_rand.argmax()]
        max_acq = y_rand.max()

        # Use LBFGSB to interpolate max
        x_seeds = random_state.uniform(bounds[:, 0], bounds[:, 1],
                                       size=(n_iter, bounds.shape[0]))
        for x_try in x_seeds:
            # Find the minimum of minus the acquisition function
            res = minimize(lambda x: -utility_func(x.reshape(1, -1), gp=gp, y_max=y_max),
                           x_try.reshape(1, -1),
                           bounds=bounds,
                           method="L-BFGS-B")

            # See if success
            if not res.success:
                continue

            # Store it if better than previous minimum(maximum).
            if max_acq is None or -res.fun[0] >= max_acq:
                x_max = res.x
                max_acq = -res.fun[0]

        # Clip output to make sure it lies within the bounds. Due to floating
        # point technicalities this is not always the case.
        return np.clip(x_max, bounds[:, 0], bounds[:, 1])


class UtilityFunction:
    def __init__(self, function, hyperparam):
        self.hyperparam = hyperparam

        if function not in ['ucb', 'ei', 'poi']:
            err = 'Utility function not implemented.' \
                  'Please choose between ucb, ei or poi'
            raise NotImplementedError(err)
        else:
            self.function = function

    def utility(self, space, gp, y_max):
        if self.function == 'ucb':
            return self._ucb(space, gp, self.hyperparam)
        if self.function == 'ei':
            return self._ei(space, gp, y_max, self.hyperparam)
        if self.function == 'poi':
            return self._poi(space, gp, y_max, self.hyperparam)

    @staticmethod
    def _ucb(space, gp, kappa):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mean, std = gp.predict(space, return_std=True)

        return mean + kappa * std

    @staticmethod
    def _ei(space, gp, y_max, xi):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mean, std = gp.predict(space, return_std=True)

        z = (mean - y_max - xi) / std
        return (mean - y_max - xi) * norm.cdf(z) + std * norm.pdf(z)

    @staticmethod
    def _poi(space, gp, y_max, xi):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mean, std = gp.predict(space, return_std=True)

        z = (mean - y_max - xi)/std
        return norm.cdf(z)
