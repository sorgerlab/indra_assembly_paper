import numpy as np
import pandas as pd
from multiprocessing import Pool
from sklearn.ensemble import RandomForestClassifier
from bioexp.curation.model_fit import ModelFit, ens_sample
from sklearn.linear_model import LogisticRegression
from bioexp.curation.belief_models import OrigBeliefStmt


class BinaryRandomForest(RandomForestClassifier):
    """Random Forest model that transforms raw counts into 0 or 1."""
    @staticmethod
    def _binarize(x_arr):
        bin_arr = x_arr.copy()
        bin_arr[bin_arr > 0] = 1
        return bin_arr

    def fit(self, x_train, y_train, *args, **kwargs):
        return super().fit(self._binarize(x_train), y_train, *args, **kwargs)

    def predict(self, x_arr, *args, **kwargs):
        return super().predict(self._binarize(x_arr), *args, **kwargs)

    def predict_proba(self, x_arr, *args, **kwargs):
        return super().predict_proba(self._binarize(x_arr), *args, **kwargs)


class LogLogisticRegression(LogisticRegression):
    """Logistic regression model that log-transforms the counts data."""
    def fit(self, x_train, y_train, *args, **kwargs):
        return super().fit(np.log(x_train+1), y_train, *args, **kwargs)

    def predict(self, x_arr, *args, **kwargs):
        return super().predict(np.log(x_arr+1), *args, **kwargs)

    def predict_proba(self, x_arr, *args, **kwargs):
        return super().predict_proba(np.log(x_arr+1), *args, **kwargs)


class BeliefModel(object):
    """Wrapper of belief models implementing sklearn classifier interface.

    reader_list : list
        List of sources.
    model_class : class or None
        One of the belief models in bioexp.curation.belief_models. If not
        provided, OrigBeliefStmt (original two-parameter Belief Model) is
        used.
    nwalkers : int
        Number of MCMC walkers.
    burn_steps : int
        Number of MCMC burn-in steps.
    sample_steps : int
        Number of MCMC sampling steps.
    """
    def __init__(self, reader_list, model_class=None, nwalkers=100,
                  burn_steps=100, sample_steps=100):
        if model_class is None:
            model_class = OrigBeliefStmt
        self.reader_list = reader_list
        self.model_class = model_class
        self.nwalkers = nwalkers
        self.burn_steps = burn_steps
        self.sample_steps = sample_steps
        self.reader_results = {}

    @staticmethod
    def df_to_num_ev(df):
        d = {}
        for _, num_ev, correct in df.itertuples():
            if num_ev not in d:
                d[num_ev] = []
            d[num_ev].append(correct)
        return d

    def fit(self, x_train, y_train, y_target=1):
        data = np.column_stack([x_train, y_train])
        self.y_ix = data.shape[1]-1
        self.y_target = y_target
        cols = self.reader_list + ['correct']
        df = pd.DataFrame(data, columns=cols)
        # Get the unique input vectors in x_train
        for reader in self.reader_list:
            r_df = df[df[reader] > 0][[reader, 'correct']]
            correct_by_num_ev = self.df_to_num_ev(r_df)
            # Convert the dataframe into a dictionary of corrects and
            # incorrects keyed by numbers of evidences
            print(reader, r_df.shape)
            model = OrigBeliefStmt()
            mf = ModelFit(model, correct_by_num_ev)
            with Pool() as pool:
                sampler = ens_sample(mf, self.nwalkers, self.burn_steps,
                                     self.sample_steps, pool=pool)
            self.reader_results[reader] = (mf, sampler)

    def predict_proba(self, x_arr):
        y_probs = np.zeros((x_arr.shape[0], 2))
        reader_errs = np.zeros((x_arr.shape[0], len(self.reader_list)))
        for ix, reader in enumerate(self.reader_list):
            x_data = x_arr[:, ix]
            mf, sampler = self.reader_results[reader]
            map_params_dict = mf.get_map_params(sampler)
            params = [map_params_dict[pname] for pname in mf.model.param_names]
            reader_errs[:, ix] = mf.model.stmt_predictions(params, x_data)
        err_probs = 1 - reader_errs
        y_probs[:, 0] = err_probs.prod(axis=1)
        y_probs[:, 1] = 1 - y_probs[:, 0]
        return y_probs

    def predict(self, x_arr, threshold=0.5):
        y_preds = np.zeros(x_arr.shape[0])
        y_probs = self.predict_proba(x_arr)
        for row_ix, pred_prob in enumerate(y_probs):
            if pred_prob[1] is np.nan:
                pred = np.nan
            else:
                pred = 0 if pred_prob[1] < threshold else 1
            y_preds[row_ix] = pred
        return y_preds
