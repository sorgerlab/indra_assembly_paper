import numpy as np
#from scipy.stats import binom, betabinom
from scipy.special import betaln, comb, beta as beta_func
from .binom_funcs import binom_pmf, binom_log_pmf, \
                         betabinom_pmf, betabinom_log_pmf


__all__ = ['BinomialEv', 'BinomialStmt', 'BetaBinomialEv', 'BetaBinomialStmt',
           'OrigBeliefEv', 'OrigBeliefStmt']

class BeliefModel(object):
    def __init__(self, param_names):
        self.param_names = param_names

    def log_prior(self, params, args):
        raise NotImplementedError()

    def log_likelihood(self, params, data, args):
        raise NotImplementedError()

    def sample_prior(self):
        raise NotImplementedError()


def bernoulli_lkl(p_by_num_ev, correct_by_num_ev):
    ll = 0
    for num_ev, num_corrects in correct_by_num_ev.items():
        p = p_by_num_ev[num_ev]
        for num_correct in num_corrects:
            if num_correct == 0:
                ll += np.log(1 - p)
            else:
                ll += np.log(p)
    return ll


# BASIC BINOMIAL MODEL ---------------------------------------------------
class Binomial(BeliefModel):
    def __init__(self):
        super(Binomial, self).__init__('p')

    def log_prior(self, params, args):
        p = params[0]
        if p < 0:
            return -np.inf
        elif p > 1:
            return - np.inf
        else:
            return 0

    def log_likelihood_ev(self, params, correct_by_num_ev, args):
        p = params[0]
        ll = 0
        for num_ev, num_corrects in correct_by_num_ev.items():
            n = num_ev
            for num_correct in num_corrects:
                k = num_correct
                ll += binom_log_pmf(k, n, p)
        return ll

    def log_likelihood_stmt(self, params, correct_by_num_ev, args):
        p = params[0]
        ll = 0
        for num_ev, num_corrects in correct_by_num_ev.items():
            n = num_ev
            for num_correct in num_corrects:
                prob_zero = binom_pmf(0, num_ev, p)
                prob_non_zero = 1 - prob_zero
                if num_correct == 0:
                    ll += np.log(prob_zero)
                else:
                    ll += np.log(prob_non_zero)
        return ll

    def sample_prior(self):
        return np.random.random(size=1)

    def ev_predictions(self, params, n):
        # Return the vector of probabilities of exactly 0 <= k <= n evidences
        # correct
        p = params[0]
        probs = []
        for k in range(0, n+1):
            probs.append(binom_pmf(k, n, p))
        return probs

    def stmt_predictions(self, params, num_evs):
        # Return the vector of probabilities correctness for statements with
        # different numbers of evidences
        p = params[0]
        probs = []
        for num_ev in num_evs:
            probs.append(1 - binom_pmf(0, num_ev, p))
        return probs


class BinomialStmt(Binomial):
    def log_likelihood(self, *args):
        return self.log_likelihood_stmt(*args)


class BinomialEv(Binomial):
    def log_likelihood(self, *args):
        return self.log_likelihood_ev(*args)


# BETA-BINOMIAL MODEL ---------------------------------------------------
class BetaBinomial(BeliefModel):
    def __init__(self):
        super(BetaBinomial, self).__init__(['Alpha', 'Beta'])

    def log_prior(self, params, args):
        # alpha and beta are positive real numbers
        alpha, beta = params
        if alpha < 0 or beta < 0:
            return -np.inf
        # FIXME: Restricting to < 1 due to NaN
        elif alpha > 1 or beta > 1:
            return -np.inf
        else:
            return 0

    def _log_lkl_k_ev(self, k, n, alpha, beta):
        b1 = betaln(k + alpha, n - k + beta)
        b2 = betaln(alpha, beta)
        nck = -betaln(1 + n - k, 1 + k) - np.log(n + 1)
        return nck + b1 - b2

    def _lkl_k_ev(self, k, n, alpha, beta):
        b1 = beta_func(k + alpha, n - k + beta)
        b2 = beta_func(alpha, beta)
        return comb(n, k) * b1 / b2

    def log_likelihood_ev(self, params, correct_by_num_ev, args):
        ll = 0
        alpha, beta = params
        for num_ev, num_corrects in correct_by_num_ev.items():
            for num_correct in num_corrects:
                #ll += self._log_lkl_k_ev(num_correct, num_ev, alpha, beta)
                ll += betabinom_log_pmf(num_correct, num_ev, alpha, beta)
        return ll

    def log_likelihood_stmt(self, params, correct_by_num_ev, args):
        ll = 0
        alpha, beta = params
        for num_ev, num_corrects in correct_by_num_ev.items():
            for num_correct in num_corrects:
                #prob_zero = self._lkl_k_ev(0, num_ev, alpha, beta)
                prob_zero = betabinom_pmf(0, num_ev, alpha, beta)
                if num_correct == 0:
                    ll += np.log(prob_zero)
                else:
                    ll += np.log(1 - prob_zero)
        return ll

    def sample_prior(self):
        # alpha and beta are positive real numbers--can scale to random reals
        # in arbitrary interval
        return np.random.random(size=2)# * 10

    def ev_predictions(self, params, n):
        # Return the vector of probabilities of exactly 0 <= k <= n evidences
        # correct
        alpha, beta = params
        probs = []
        for k in range(0, n+1):
            #probs.append(self._lkl_k_ev(k, n, alpha, beta))
            probs.append(betabinom_pmf(k, n, alpha, beta))
        return probs

    def stmt_predictions(self, params, num_evs):
        # Return the vector of probabilities correctness for statements with
        # different numbers of evidences from 1 to max_n
        alpha, beta = params
        probs = []
        for num_ev in num_evs:
            #prob_zero = self._lkl_k_ev(0, num_ev, alpha, beta)
            prob_zero = betabinom_pmf(0, num_ev, alpha, beta)
            probs.append(1 - prob_zero)
        return probs


class BetaBinomialStmt(BetaBinomial):
    def log_likelihood(self, *args):
        return self.log_likelihood_stmt(*args)


class BetaBinomialEv(BetaBinomial):
    def log_likelihood(self, *args):
        return self.log_likelihood_ev(*args)


# ORIGINAL BELIEF MODEL -------------------------------------
class OrigBelief(BeliefModel):
    def __init__(self):
        super(OrigBelief, self).__init__(['Rand', 'Syst'])

    def log_prior(self, params, args):
        pr, ps = params
        # Uniform prior over [0, 1]
        if pr < 0 or ps < 0 or pr > 1 or ps > 1:
            return -np.inf
        else:
            return 0

    @staticmethod
    def belief(num_ev, pr, ps):
        b = 1 - (ps + (1-ps) * (pr ** num_ev))
        return b

    def log_likelihood_ev(self, params, correct_by_num_ev, args):
        pr, ps = params
        ll = 0
        for num_ev, num_corrects in correct_by_num_ev.items():
            for num_correct in num_corrects:
                if num_correct == 0:
                    ll += np.log(ps + (1-ps) * binom_pmf(0, num_ev, 1-pr))
                else:
                    ll += np.log((1-ps) * binom_pmf(num_correct, num_ev, 1-pr))
        return ll

    def log_likelihood_stmt(self, params, correct_by_num_ev, args):
        pr, ps = params
        ll = 0
        for num_ev, corrects in correct_by_num_ev.items():
            ll += sum((np.log(self.belief(num_ev, pr, ps)) if c else
                       np.log(1-self.belief(num_ev, pr, ps)))
                       for c in corrects)
        return ll

    def sample_prior(self):
        # pr and ps are between [0, 1)
        return np.random.random(size=2)

    def ev_predictions(self, params, n):
        """Return log likelihood of belief model parameters given data."""
        pr, ps = params
        probs = []
        for k in range(0, n+1):
            if k == 0:
                ll = ps + (1-ps) * binom_pmf(0, n, 1-pr)
            else:
                ll = (1-ps) * binom_pmf(k, n, 1-pr)
            probs.append(ll)
        return probs

    def stmt_predictions(self, params, num_evs):
        # Return the vector of probabilities correctness for statements with
        # different numbers of evidences from 1 to max_n
        pr, ps = params
        probs = []
        for num_ev in num_evs:
            probs.append(self.belief(num_ev, pr, ps))
        return probs

# ORIGINAL BELIEF MODEL by evidence -------------------------------------

class OrigBeliefStmt(OrigBelief):
    def log_likelihood(self, *args):
        return self.log_likelihood_stmt(*args)

class OrigBeliefEv(OrigBelief):
    def log_likelihood(self, *args):
        return self.log_likelihood_ev(*args)
