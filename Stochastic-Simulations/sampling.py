import math
import random as rnd
import numpy as np
import matplotlib.pyplot as plt


class ProbabilityDistribution:
    """
    A base class for representing a probability distribution, either continuous or discrete.

    Attributes:
        func (callable): A probability density/mass function.
        interval (tuple): The interval of the distribution domain (a, b).
    """

    def __init__(self, func: callable, interval: tuple):
        self.func = func
        self.interval = interval

    def get_random_x(self):
        """Generate a random x in the domain interval of the distribution."""
        a, b = self.interval
        return rnd.uniform(a, b)

    def get_random_y(self, x, start_codomain: int = 0):
        """Generate a random y for a given x, with the y values constrained by the function."""
        return rnd.uniform(start_codomain, self.f(x))

    def get_random_y_in_interval(self, codomain_interval: tuple):
        """Generate a random y in the given codomain interval."""
        a, b = codomain_interval
        return rnd.uniform(a, b)

    def f(self, x):
        """Evaluate the probability function at point x."""
        return self.func(x)


class ProbabilityDistributionContinuous(ProbabilityDistribution):
    """
    A class for representing continuous probability distributions.
    Inherits from the ProbabilityDistribution class.
    """

    def get_random_x(self):
        """Generate a continuous random x in the domain interval."""
        a, b = self.interval
        return rnd.uniform(a, b)

    def get_random_y(self, x, start_codomain: int = 0):
        """Generate a continuous random y value for a given x."""
        return rnd.uniform(start_codomain, self.f(x))

    def get_random_y_in_interval(self, codomain_interval: tuple):
        """Generate a continuous random y in the given codomain interval."""
        a, b = codomain_interval
        return rnd.uniform(a, b)


class ProbabilityDistributionDiscrete(ProbabilityDistribution):
    """
    A class for representing discrete probability distributions.
    Inherits from the ProbabilityDistribution class.
    """

    def get_random_x(self):
        """Generate a discrete random x in the domain interval."""
        a, b = self.interval
        return rnd.randrange(a, b)

    def get_random_y(self, x, start_codomain: int = 0):
        """Generate a discrete random y value for a given x."""
        return rnd.uniform(start_codomain, self.f(x))

    def get_random_y_in_interval(self, codomain_interval: tuple):
        """Generate a discrete random y in the given codomain interval."""
        a, b = codomain_interval
        return rnd.uniform(a, b)


def rejection_sampling(g_proxy: ProbabilityDistribution, f: callable) -> tuple:
    while True:
        x_sampling = g_proxy.get_random_x()
        f_xo = f(x_sampling)
        g_xo = g_proxy.f(x_sampling)

        y_sampling = g_proxy.get_random_y_in_interval((0, g_xo))

        if y_sampling <= f_xo:
            return x_sampling, y_sampling


def get_random_from_poisson(_lambda: float):
    if _lambda < 10:
        max = (_lambda ** _lambda * np.exp(-_lambda)) / math.factorial(int(_lambda))
        # approximation to _lambda * 3 as last possible X value
        uniform_dist = ProbabilityDistributionDiscrete(lambda x: max, (0, int(np.ceil(_lambda))))
        return rejection_sampling(uniform_dist, lambda x: (_lambda ** x * np.exp(-_lambda)) / math.factorial(int(x)))
    else:
        # approximation with gaussian distribution
        return get_random_from_gaussian(_lambda,_lambda)


def get_random_from_gaussian(mu, std, interval: tuple = (-100, 100)):
    def normal_distribution(mu, std, x):
        return 1 / np.sqrt(2 * np.pi) * std * np.exp(-1 / 2 * ((x - mu) / std) ** 2)

    max = normal_distribution(mu, std, mu)

    uniform_dist = ProbabilityDistributionDiscrete(lambda x: max, interval)

    return rejection_sampling(uniform_dist, lambda x: normal_distribution(mu, std, x))


