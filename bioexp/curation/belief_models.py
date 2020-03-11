import numpy as np
from scipy.stats import binom, betabinom
from scipy.special import betaln, comb, beta as beta_func


class BeliefModel(object):
    def __init__(self, param_names):
        self.param_names = param_names

    def log_prior(self, params, args):
        raise NotImplementedError()

    def log_likelihood(self, params, data, args):
        raise NotImplementedError()

    def sample_prior(self):
        raise NotImplementedError()



# BASIC BINOMIAL MODEL ---------------------------------------------------
class BinomialModelEv(BeliefModel):
    def __init__(self):
        super(BinomialModelEv, self).__init__('p')

    def log_prior(self, params, args):
        p = params[0]
        if p < 0:
            return -np.inf
        elif p > 1:
            return - np.inf
        else:
            return 0

    @staticmethod
    def _log_lkl_k_ev(k, n, p):
        """Binomial PMF: (n choose k) * p^k * (1 - p)^(n-k)"""
        nck = -betaln(1 + n - k, 1 + k) - np.log(n + 1)
        pk = k * np.log(p)
        qnmk = (n - k) * np.log(1 - p)
        return nck + pk + qnmk

    @staticmethod
    def _lkl_k_ev(k, n, p):
        return comb(n, k) * (p ** k) * ((1 - p) ** (n - k))

    def log_likelihood(self, params, correct_by_num_ev, args):
        p = params[0]
        ll = 0
        for num_ev, num_corrects in correct_by_num_ev.items():
            n = num_ev
            for num_correct in num_corrects:
                k = num_correct
                #ll += np.log(binom.pmf(k, n, p))
                ll += self._log_lkl_k_ev(k, n, p)
        return ll

    def sample_prior(self):
        return np.random.random(size=1)

    def ev_predictions(self, params, n):
        # Return the vector of probabilities of exactly 0 <= k <= n evidences
        # correct
        p = params[0]
        probs = []
        for k in range(0, n+1):
            #probs.append(binom.pmf(k, n, p))
            probs.append(self._lkl_k_ev(k, n, p))
            #np.exp(self._log_lkl_k_ev(k, n, p)))
        return probs

    def stmt_predictions(self, params, num_evs):
        # Return the vector of probabilities correctness for statements with
        # different numbers of evidences
        p = params[0]
        probs = []
        for num_ev in num_evs:
            #probs.append(binom.sf(0, num_ev, p))
            probs.append(1 - self._lkl_k_ev(0, num_ev, p))
        return probs


# BINOMIAL by statement ---------------------------------------------
class BinomialModelStmt(BinomialModelEv):
    def log_likelihood(self, params, correct_by_num_ev, args):
        p = params[0]
        ll = 0
        for num_ev, num_corrects in correct_by_num_ev.items():
            n = num_ev
            for num_correct in num_corrects:
                prob_zero = self._lkl_k_ev(0, num_ev, p)
                prob_non_zero = 1 - prob_zero
                if num_correct == 0:
                    ll += np.log(prob_zero)
                else:
                    ll += np.log(prob_non_zero)
        return ll


# BETA-BINOMIAL MODEL ---------------------------------------------------
class BetaBinomialModelEv(BeliefModel):
    def __init__(self):
        super(BetaBinomialModelEv, self).__init__(['Alpha', 'Beta'])

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
        """Log version of the likelihood."""
        b1 = betaln(k + alpha, n - k + beta)
        b2 = betaln(alpha, beta)
        nck = -betaln(1 + n - k, 1 + k) - np.log(n + 1)
        return nck + b1 - b2

    def _lkl_k_ev(self, k, n, alpha, beta):
        """Likelihood."""
        b1 = beta_func(k + alpha, n - k + beta)
        b2 = beta_func(alpha, beta)
        return comb(n, k) * b1 / b2

    def log_likelihood(self, params, correct_by_num_ev, args):
        ll = 0
        alpha, beta = params
        for num_ev, num_corrects in correct_by_num_ev.items():
            for num_correct in num_corrects:
                ll += self._log_lkl_k_ev(num_correct, num_ev, alpha, beta)
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
            probs.append(self._lkl_k_ev(k, n, alpha, beta))
        return probs

    """
    @staticmethod
    def _cdf(k, n, alpha, beta):
        # n choose k
        if k < 0:
            return 0
        elif k >= n:
            return 1
        else:
            return (comb(n, k) * beta_func(k + alpha, n - k + beta) /
                                                    beta_func(alpha, beta))
    """

    def stmt_predictions(self, params, num_evs):
        # Return the vector of probabilities correctness for statements with
        # different numbers of evidences from 1 to max_n
        alpha, beta = params
        probs = []
        for num_ev in num_evs:
            prob_zero = self._lkl_k_ev(0, num_ev, alpha, beta)
            probs.append(1 - prob_zero)
        return probs


# BETA-BINOMIAL by statement ---------------------------------------------
class BetaBinomialModelStmt(BetaBinomialModelEv):
    def log_likelihood(self, params, correct_by_num_ev, args):
        ll = 0
        alpha, beta = params
        for num_ev, num_corrects in correct_by_num_ev.items():
            for num_correct in num_corrects:
                prob_zero = self._lkl_k_ev(0, num_ev, alpha, beta)
                if prob_zero == 1:
                    import ipdb; ipdb.set_trace()
                ll += np.log(1 - prob_zero)
        if ll == np.nan:
            import ipdb; ipdb.set_trace()
        return ll



# ORIGINAL BELIEF MODEL by statement -------------------------------------
class OrigBeliefModelStmt(BeliefModel):
    def __init__(self):
        super(OrigBeliefModelStmt, self).__init__(['Rand', 'Syst'])

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

    def log_likelihood(self, params, correct_by_num_ev, args):
        """Return log likelihood of belief model parameters given data."""
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
                ll = ps + (1-ps)*binom.pmf(0, n=n, p=1-pr)
            else:
                ll = (1-ps) * binom.pmf(k, n=n, p=1-pr)
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

class OrigBeliefModelEv(OrigBeliefModelStmt):
    @staticmethod
    def belief(num_ev, pr, ps):
        b = 1 - (ps + (1-ps) * (pr ** num_ev))
        return b

    def binom_pmf(k, n, p):
        return comb(n, k) * p**k * (1-p)**(n-k)

    def log_likelihood(self, params, correct_by_num_ev, args):
        """Return log likelihood of belief model parameters given data."""
        pr, ps = params
        ll = 0
        for num_ev, num_corrects in correct_by_num_ev.items():
            for num_correct in num_corrects:
                if num_correct == 0:
                    ll += np.log(ps + (1-ps) *
                                        self.binom_pmf(0, n, p=1-pr))
                else:
                    ll += np.log((1-ps) *
                               self.binom_pmf(num_correct, num_ev, 1-pr))
        return ll
