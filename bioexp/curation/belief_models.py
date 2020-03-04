import numpy as np
from scipy.stats import binom
from scipy.special import betaln
from bioexp.curation.model_fit import BeliefModel

# BASIC BINOMIAL MODEL ---------------------------------------------------
class BinomialModelEv(BeliefModel):
    def log_prior(self, params, args):
        p = params[0]
        if p < 0:
            return -np.inf
        elif p > 1:
            return - np.inf
        else:
            return 0

    def log_likelihood(self, params, correct_by_num_ev, args):
        p = params[0]
        ll = 0
        for num_ev, num_corrects in correct_by_num_ev.items():
            n = num_ev
            for num_correct in num_corrects:
                k = num_correct
                ll += np.log(binom.pmf(k, n, p))
        return ll

    def sample_prior(self):
        return np.random.random()


    def ev_predictions(self, params, n):
        # Return the vector of probabilities of exactly 0 <= k <= n evidences
        # correct
        p = params[0]
        probs = []
        for k in range(0, n+1):
            probs.append(binom.pmf(k, n, p))
        return probs

    def stmt_predictions(self, params, num_evs):
        # Return the vector of probabilities correctness for statements with
        # different numbers of evidences
        p = params[0]
        probs = []
        for num_ev in num_evs:
            probs.append(binom.sf(0, num_ev, p))
        return probs

# BETA-BINOMIAL MODEL ---------------------------------------------------
class BetaBinomialModelEv(BeliefModel):
    def log_prior(self, params, args):
        # alpha and beta are positive real numbers
        alpha, beta = params
        if alpha < 0 or beta < 0:
            return -np.inf
        else:
            return 0

    def _log_lkl_k_ev(self,k, n, alpha, beta):
        k = num_correct
        b1 = betaln(k + alpha, n - k + beta)
        b2 = betaln(alpha, beta)
        nck = -betaln(1 + n - k, 1 + k) - np.log(n + 1)
        return nck + b1 - b2

    def log_likelihood(self, params, correct_by_num_ev, args):
        ll = 0
        alpha, beta = params
        for num_ev, num_corrects in correct_by_num_ev.items():
            for num_correct in num_corrects:
                ll += self._log_lkl_k_ev(num_correct, num_ev, alpha, beta)
        return ll

    def sample_prior(self):
        # alpha and beta are positive real numbers--scale to random reals
        # in arbitrary interval [0, 10)
        return np.random.random(size=2) * 10

    def ev_predictions(self, params, n):
        # Return the vector of probabilities of exactly 0 <= k <= n evidences
        # correct
        alpha, beta = params
        probs = []
        for k in range(0, n+1):
            probs.append(self._log_lkl_k_ev(k, n, alpha, beta))
        return probs

    def stmt_predictions(self, params, num_evs):
        # Return the vector of probabilities correctness for statements with
        # different numbers of evidences from 1 to max_n
        p = params[0]
        probs = []
        for num_ev in num_evs:
            probs.append(binom.sf(1, num_ev, p))
        return probs


