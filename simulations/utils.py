import numpy as np
import pickle

TX_STRAND_DICT = {'strand1':range(0, 10), 'strand2':range(10, 19)}

def beta_rv(alpha, beta, **kwargs):
    return np.random.beta(a=alpha, b=beta, **kwargs)

def binomial_rv(n_trials, prob, **kwargs):
    return np.random.binomial(n=n_trials, p=prob, **kwargs)

def geometric_rv(n_sucesses, prob):
    return np.random.geometric(p=prob, size=n_sucesses)

def neg_binomial_rv(success_number, prob):
    ## formulation is number of failures till that success
    return np.random.negative_binomial(n=success_number, p=prob)

def gamma_rv(event_number, rate):
    return np.random.gamma(shape=event_number, scale=1/rate)

def poisson_rv(rate, exposure):
    return np.random.poisson(lam=rate*exposure)

def normal_rv(mean, std_dev, **kwargs):
    return np.random.normal(loc=mean, scale=std_dev, **kwargs)

