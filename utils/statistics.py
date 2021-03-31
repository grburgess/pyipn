import scipy.stats as stats
from scipy.stats import norm
from scipy.optimize import fsolve


def p2sigma(p):
    def f(x):
        return (norm.cdf(x[0]) - norm.cdf(-x[0])) - p

    return fsolve(f, [1])[0]


def sigma2prob(sigma):
    return stats.norm.cdf(sigma) - stats.norm.cdf(-sigma)
