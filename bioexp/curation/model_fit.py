import numpy as np
from matplotlib import pyplot as plt
import emcee
from texttable import Texttable
import corner

def posterior(position, mf):
    """A generic log posterior function."""
    pr = prior(position, mf)
    if pr == -np.inf:
        return -np.inf
    else:
        return pr + likelihood(position, mf)

def prior(position, mf):
    """A generic prior function."""
    return mf.model.log_prior(position, mf)

def likelihood(position, mf):
    return mf.model.log_likelihood(position, mf.data, None)



class ModelFit(object):
    def __init__(self, model, data):
        self.model = model
        self.data = data
        self.data_stmt = {}
        # Convert at the evidence level to data at the stmt level and store
        for num_ev in data.keys():
            stmt_corrects = list(np.array(np.array(data[num_ev]) >= 1,
                                          dtype=int))
            self.data_stmt[num_ev] = stmt_corrects

    def plot_ev_fit(self, sampler):
        plt.figure()
        map_ix = np.argmax(sampler.flatlnprobability)
        map_p = sampler.flatchain[map_ix]
        for n in range(1, 11):
            # First, plot the data
            plt.subplot(3, 4, n)
            plt.hist(self.data[n], bins=range(0,n+2), density=True)
            plt.title(n)
            # Then plot the model predictions
            ev_lks = self.model.ev_predictions(map_p, n)
            bin_centers = 0.5 + np.array(range(0, n+1))
            print(bin_centers)
            plt.plot(bin_centers, ev_lks, color='r', marker='.')
        plt.tight_layout()
        plt.show()

    def plot_stmt_fit(self, sampler):
        # Calculate the mean of correctness by number of evidence
        num_evs = sorted(self.data.keys())
        means = [np.mean(self.data_stmt[n]) for n in num_evs]
        # Stderr of proportion is sqrt(pq/n)
        std = [2*np.sqrt((np.mean(self.data_stmt[n]) *
                          (1 - np.mean(self.data_stmt[n]))) /
                           len(self.data_stmt[n]))
               for n in num_evs]
        #beliefs = [belief(n, opt_r, opt_s) for n in num_evs]

        # Print table of results before plotting
        table = Texttable()
        table_data = [['Num Evs', 'Count', 'Num Correct', 'Pct', 'Std']]
        for i, num_ev in enumerate(num_evs):
            table_row = [num_ev, len(self.data_stmt[num_ev]),
                         sum(self.data_stmt[num_ev]), means[i], std[i]]
            table_data.append(table_row)
        table.add_rows(table_data)
        print(table.draw())

        # Plot the data
        plt.figure()
        plt.errorbar(num_evs, means, yerr=std, fmt='bo-', ls='none',
                     label='Empirical mean correctness')
        # Plot the MAP predictions
        map_ix = np.argmax(sampler.flatlnprobability)
        map_p = sampler.flatchain[map_ix]
        stmt_probs = self.model.stmt_predictions(map_p, num_evs)
        plt.plot(num_evs, stmt_probs, 'ro-', label='MAP belief')
        # Legends, labels, etc.
        plt.ylim(0, 1)
        plt.grid(True)
        plt.xticks(num_evs)
        plt.xlabel('Number of evidence per INDRA Statement')
        plt.legend(loc='lower right')
        plt.show()

    def plot_corner(self, sampler):
        # Plot the posterior parameter distribution
        plt.figure()
        corner.corner(sampler.flatchain, labels=self.model.param_names)
        plt.show()
        """
        # Plot a few representative belief curves from the posterior
        num_evs = sorted(ev_correct_by_num_ev.keys())
        for pr, ps in samples[:100]:
            beliefs = [belief(n, pr, ps) for n in num_evs]
            plt.plot(num_evs, beliefs, 'g-', alpha=0.1)
        """

def ens_sample(mf, nwalkers, burn_steps, sample_steps, threads=1,
               pos=None, random_state=None, pool=None):
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
                p0[walk_ix, :] = mf.model.sample_prior()
    else:
        p0 = pos

    # Create the sampler object
    sampler = emcee.EnsembleSampler(nwalkers, ndim, posterior,
                                         args=[mf],
                                         threads=threads, pool=pool)
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

