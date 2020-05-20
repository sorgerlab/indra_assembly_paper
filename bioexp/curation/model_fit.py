import emcee
import corner
import logging
import numpy as np
from texttable import Texttable
from matplotlib import pyplot as plt


logger = logging.getLogger('model_fit')


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

    def data_table(self):
        # Print table of results before plotting
        table = Texttable()
        table_data = [['Num Evs', 'Count', 'Num Correct', 'Pct', 'Std']]
        for i, num_ev in enumerate(num_evs):
            table_row = [num_ev, len(self.data_stmt[num_ev]),
                         sum(self.data_stmt[num_ev]), means[i], std[i]]
            table_data.append(table_row)
        table.add_rows(table_data)
        print(table.draw())

    def stmt_err(self, sampler, weights=None):
        map_ix = np.argmax(sampler.flatlnprobability)
        map_p = sampler.flatchain[map_ix]
        stmt_preds = self.model.stmt_predictions(map_p, self.data.keys())
        ll = 0
        for i, (num_ev, num_corrects) in enumerate(self.data.items()):
            p = stmt_preds[i]
            ll_n = 0
            for num_correct in num_corrects:
                if num_correct == 0:
                    ll_n += np.log(1 - p)
                else:
                    ll_n += np.log(p)
            if weights:
                ll += weights[num_ev] * ll_n * len(self.data)
            else:
                ll += ll_n
        return -ll

    def plot_ev_fit(self, sampler, title):
        fig = plt.figure()
        map_ix = np.argmax(sampler.flatlnprobability)
        map_p = sampler.flatchain[map_ix]
        for n in range(1, max(self.data.keys())+1):
            # First, plot the data
            plt.subplot(3, 4, n)
            if n not in self.data:
                continue
            plt.hist(self.data[n], bins=range(0, n+2), density=True)
            plt.title(n)
            # Then plot the model predictions
            ev_lks = self.model.ev_predictions(map_p, n)
            bin_centers = 0.5 + np.array(range(0, n+1))
            plt.plot(bin_centers, ev_lks, color='r', marker='.')
        fig.suptitle(title)
        plt.tight_layout()
        plt.show()

    def plot_stmt_fit(self, sampler, title):
        # Calculate the mean of correctness by number of evidence
        num_evs = sorted(self.data.keys())
        means = [np.mean(self.data_stmt[n]) for n in num_evs]
        # Stderr of proportion is sqrt(pq/n)
        std = [2*np.sqrt((np.mean(self.data_stmt[n]) *
                          (1 - np.mean(self.data_stmt[n]))) /
                           len(self.data_stmt[n]))
               for n in num_evs]

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
        plt.title(title)
        plt.show()

    def plot_corner(self, sampler):
        # Plot the posterior parameter distribution
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

    logger.info("Burn in sampling...")
    pos, prob, state = sampler.run_mcmc(p0, burn_steps, store=False)
    sampler.reset()

    logger.info("Main sampling...")
    sampler.run_mcmc(pos, sample_steps)

    logger.info("Done sampling.")
    return sampler

