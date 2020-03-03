import numpy as np
import emcee

def posterior(position, mf):
    """A generic log posterior function."""
    pr = prior(position, mf)
    if pr == -np.inf:
        return -np.inf
    else:
        return pr + likelihood(position, mf)

def prior(position, mf):
    """A generic prior function."""
    return mf.model.log_prior(position)

def likelihood(position, mf):
    return mf.model.log_likelihood(position, mf.data)


class BeliefModel(object):
    def __init__(self, param_names):
        self.param_names = param_names

    def log_prior(self, params):
        raise NotImplementedError()

    def log_likelihood(self, params, data):
        raise NotImplementedError()

    def sample_prior(self):
        raise NotImplementedError()


class ModelFit(object):
    def __init__(self, model, data):
        self.model = model
        self.data = data

    def plot_map_fit(self, sampler):
        # FIXME: Plotting function specific to belief data here
        plt.figure()
        plt.plot(range(len(data)), self.data)
        ml_ix = np.argmax(sampler.flatlnprobability)
        ml_p = sampler.flatchain[ml_ix]
        plt.plot(range(len(data)), posterior(

def ens_sample(mf, nwalkers, burn_steps, sample_steps, threads=1,
               pos=None, random_state=None):
    """Samples from the posterior function using emcee.EnsembleSampler.

    The EnsembleSampler containing the chain is stored in gf.sampler.

    Note that parameters are log10-transformed during fitting, so the
    parameter values returned from the walk must be exponentiated to get
    them back to untransformed values (e.g., 10 ** gf.sampler.flatchain)

    Parameters
    ----------
    mf : ModelFit
        GlobalFit object containing the timepoints, data, builder object,
        Solver, etc.
    nwalkers : int
        Number of walkers to use in the emcee sampler.
    burn_steps : int
        Number of burn-in steps.
    sample_steps : int
        Number of sampling steps.
    threads : int
        Number of threads to use for parallelization. Default is 1.
    pos : numpy.array
        Matrix of initial positions for the chain. If None (default) random
        positions are chosen from the prior. Assigning a position allows
        previously run chains to be extended.
    random_state : random state for Mersenne Twister PRNG
        The random state to use to initialize the sampler's pseudo-random
        number generator. Can be used to continue runs from previous ones.
    """
    # Initialize the parameter array with initial values (in log10 units)
    # Number of parameters to estimate
    ndim = len(mf.model.param_names)
    # Initialize the walkers with starting positions drawn from the priors
    # Note that the priors are in log10 scale already, so they don't
    # need to be transformed here
    if pos is None:
        p0 = np.zeros((nwalkers, ndim))
        for walk_ix in range(nwalkers):
            for p_ix in range(ndim):
                p0[walk_ix, p_ix] = mf.model.sample_prior()
    else:
        p0 = pos

    # Create the sampler object
    sampler = emcee.EnsembleSampler(nwalkers, ndim, posterior,
                                         args=[mf],
                                         threads=threads)
    if random_state is not None:
        sampler.random_state = random_state

    print("Burn in sampling...")
    pos, prob, state = sampler.run_mcmc(p0, burn_steps, storechain=False)
    sampler.reset()

    print("Main sampling...")
    sampler.run_mcmc(pos, sample_steps)

    print("Done sampling.")
    return sampler

if __name__ == '__main__':
    data = np.random.randint(2, size=1000)

    class BernoulliModel(BeliefModel):
        def log_prior(self, params):
            p = params[0]
            if p < 0:
                return -np.inf
            elif p > 1:
                return - np.inf
            else:
                return 0

        def log_likelihood(self, params, data):
            p = params[0]
            ll = sum(np.log(p) if d == 0 else np.log(1 - p)
                     for d in data)
            return ll

        def sample_prior(self):
            return np.random.random()


    bm = BernoulliModel(['p'])
    mf = ModelFit(bm, data)

    nwalkers, burn_steps, sample_steps = (100, 1000, 100)
    sampler = ens_sample(mf, nwalkers, burn_steps, sample_steps)

