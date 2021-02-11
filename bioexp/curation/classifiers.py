import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression

class BinaryRandomForest(RandomForestClassifier):
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
    def fit(self, x_train, y_train, *args, **kwargs):
        return super().fit(np.log(x_train+1), y_train, *args, **kwargs)

    def predict(self, x_arr, *args, **kwargs):
        return super().predict(np.log(x_arr+1), *args, **kwargs)

    def predict_proba(self, x_arr, *args, **kwargs):
        return super().predict_proba(np.log(x_arr+1), *args, **kwargs)
