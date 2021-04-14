import numpy as np
import pandas as pd
import scipy.optimize
# import sklearn.linear_model


class LinearModel(object):
    def __init__(self, **kwds):
        self.coef = None

    def __call__(self, Xnew):
        if self.coef is None:
            return np.zeros(len(Xnew))
        else:
            return np.dot(Xnew, self.coef)

    def predict(self, Xnew):
        return self(Xnew)


class NNLS(LinearModel):
    def fit(self, X, Y):
        coef, residual_sqrt = scipy.optimize.nnls(X.values, Y.values)
        self.coef = pd.Series(coef, index=X.columns)
        self.l2 = residual_sqrt**2
        self.fitted = pd.Series(self.predict(X), index=Y.index)
        self.res = Y - self.fitted


class DeconvolvedUnderfit(LinearModel):
    def fit(self, X, Y, lam=1.0):
        X_nonzero = X.loc[X.index >= 0]
        Y_nonzero = Y.loc[Y.index >= 0]
        N, D = X_nonzero.shape
        fitted_prob = X_nonzero.sum(axis=0)
        non_fitted_prob = X.iloc[X.index==-1].sum(axis=0)
        lam = min((fitted_prob / non_fitted_prob).min(), lam)
        c = np.vstack( lam*np.ones(D) - (lam+1)*fitted_prob )
        self.res = scipy.optimize.linprog(c=c, A_ub=X_nonzero.values, b_ub=Y_nonzero.values)
        self.fun = Y_nonzero.sum() - self.res['fun']
        self.coef = pd.Series(self.res['x'], index=X_nonzero.columns)
        self.fitted = pd.Series(self.predict(X_nonzero), index=Y_nonzero.index)
        self.res = Y - self.fitted


class SKLEARN_WRAPPER(LinearModel):
    def fit(self, X, Y):
        self.model.fit(X,Y)
        self.coef = pd.Series(self.model.coef_, index=X.columns)
        self.fitted = pd.Series(self.model.predict(X), index=Y.index)
        self.res = Y - self.fitted
        self.intercept = self.model.intercept_


# class sklearnNNLS(SKLEARN_WRAPPER):
#     def __init__(self,
#                  fit_intercept=False,
#                  normalize=False,
#                  copy_X=False,
#                  n_jobs=-1,
#                  positive=True,
#                  **kwds):
#         self.model = sklearn.linear_model.LinearRegression(fit_intercept=fit_intercept,
#                                                            normalize=normalize,
#                                                            copy_X=copy_X,
#                                                            n_jobs=n_jobs,
#                                                            positive=positive,
#                                                            **kwds)
#         self.intercept = 0.0


# class sklearnLasso(SKLEARN_WRAPPER):
#     def __init__(self,
#                  alpha=1,
#                  fit_intercept=False,
#                  normalize=False,
#                  copy_X=False,
#                  positive=True,
#                  **kwds):
#         self.model = sklearn.linear_model.Lasso(fit_intercept=fit_intercept,
#                                                 normalize=normalize,
#                                                 copy_X=copy_X,
#                                                 positive=positive,
#                                                 **kwds)
#         self.intercept = 0.0




# class CoefConstrainedLinearRegression(LinearModel):
#     def fit(self, X, Y, upper_estimates, lam=10.0):
#         import cvxopt
#         cvxopt.solvers.options['show_progress'] = False
#         N,D = X.shape
#         P = X.values.T.dot(X)
#         q = np.vstack( 2*( -X.values.T.dot(Y) + lam*np.ones(D)) )
#         Gmat = np.block([[-np.identity(D) ],
#                          [ np.identity(D) ]])
#         beta_max = upper_estimates.loc[X.columns]
#         h = np.block([[ np.vstack(np.zeros(D))      ],
#                       [ np.vstack(beta_max.values) ]])
#         cmat = cvxopt.matrix
#         res = cvxopt.solvers.qp(P=cmat(P), q=cmat(q), G=cmat(Gmat), h=cmat(h))
#         self.coef = pd.Series(np.array(res['x']).flatten(), index=X.columns)
#         self.fitted = pd.Series(X.dot(self.coef), index=Y.index)
#         self.res = Y - self.fitted



def fit(Model):
    def fit_model(X, Y, fit_kwds={}, **kwds):
        model = Model(**kwds)
        model.fit(X,Y,**fit_kwds)
        return model
    return fit_model


fit_NNLS = fit(NNLS)
# fit_sklearnNNLS = fit(sklearnNNLS)
# fit_sklearnLasso = fit(sklearnLasso)
# fit_CCLR = fit(CoefConstrainedLinearRegression)
fit_DeconvolvedUnderfit = fit(DeconvolvedUnderfit)