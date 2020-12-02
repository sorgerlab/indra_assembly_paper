import numpy as np
from scipy.special import betaln, comb, loggamma, beta as beta_func


def log_binom(n, k):
    return loggamma(n+1) - loggamma(k+1) - loggamma(n-k+1)


def logQ(z, n):
    return (n-0.5)*np.log1p(n/z) + z*(np.log1p(n/z) - n/z)


def loggamma_star(a):
    if a == 0:
        output = np.float('inf')
    elif a > 1e9:
        output = 0.0
    else:
        output = (loggamma(a) + a - 0.5*np.log(2*np.pi) -
                  (a - 0.5)*np.log(a))
    return output


def loggamma_ratio(z, n):
    return n*np.log(z) + loggamma_star(z+n) - loggamma_star(z) + logQ(z, n)


def betabinom_log_pmf(k, n, a, b):
    log_out = (log_binom(n, k) + loggamma_ratio(a, k) + loggamma_ratio(b, n-k) -
               loggamma_ratio(a+b, n))
    return log_out


def betabinom_pmf(k, n, a, b):
    return np.exp(betabinom_log_pmf(k, n, a, b))


def binom_pmf(k, n, p):
    return comb(n, k) * (p ** k) * ((1 - p) ** (n - k))


def binom_log_pmf(k, n, p):
    nck = -betaln(1 + n - k, 1 + k) - np.log(n + 1)
    pk = k * np.log(p)
    qnmk = (n - k) * np.log(1 - p)
    return nck + pk + qnmk
