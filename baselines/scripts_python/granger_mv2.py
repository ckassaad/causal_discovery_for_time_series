import numpy as np
import pandas as pd
from statsmodels.tsa.api import VAR
from scipy.stats import f
from sklearn.preprocessing import StandardScaler


class Granger:
    def __init__(self, x, p=22, scale=True):
        """
        :param x: multivariate time series in the form of a pandas dataframe
        :param p: time stamp for prediction
        :param scale: if True normalize data
        """
        self.p = p
        self.d = x.shape[1]
        self.names = list(x.columns.values)
        self.pa = {self.names[i]: [self.names[i]] for i in range(len(self.names))}

        # scaling data
        if scale:
            scaler = StandardScaler()
            x_scaled = scaler.fit_transform(x.values)
            self.X = pd.DataFrame(x_scaled, columns=self.names)
        else:
            self.X = pd.DataFrame(x, columns=self.names)

    def predict(self, model, x):
        x_hat = pd.DataFrame(columns=x.columns.values)
        for t in range(x.shape[0] - self.p):
            temp = pd.DataFrame(model.forecast(x.values[t:(t + self.p)], 1), columns=list(x.columns.values))
            x_hat = x_hat.append(temp, ignore_index=True)
        return x_hat

    def f_test(self, var1, var2, m):
        if var1 > var2:
            f_ = np.divide(var1, var2)
        else:
            f_ = np.divide(var2, var1)
        p_value = 1 - f.cdf(f_, m - 1, m - 1)
        return p_value

    def fit(self, alpha=0.05):
        """
        :param alpha: threshold of F-test
        :return: granger causality denpendencies
        """
        model_full = VAR(self.X)
        model_full_fit = model_full.fit(maxlags=self.p, ic='aic')

        # make prediction
        x_hat = self.predict(model_full_fit, self.X)

        # compute error
        err_full =np.subtract(x_hat.values, self.X.values[self.p:])
        var_full = list(np.var(err_full, axis=0))

        for j in range(self.d):
            x_temp = self.X.drop(columns=[self.names[j]])
            model_rest = VAR(x_temp)
            model_rest_fit = model_rest.fit(maxlags=self.p, ic='aic')

            # make prediction
            x_hat = self.predict(model_rest_fit, x_temp)

            # compute error
            err_rest = np.subtract(x_hat.values, x_temp.values[self.p:])
            var_rest = list(np.var(err_rest, axis=0))

            # F test (extremely sensitive to non-normality of X and Y)
            var_full_rest = var_full.copy()
            del var_full_rest[j]
            m = x_hat.shape[0]

            for i in range(len(x_hat.columns.values)):
                # Start Test using F-test
                p_value = self.f_test(var_rest[i], var_full_rest[i], m)
                if p_value < alpha:
                    self.pa[x_hat.columns.values[i]].append(self.names[j])

        res_df = pd.DataFrame(np.ones([self.d, self.d]), columns=self.names, index=self.names)
        for e in self.pa.keys():
            for c in self.pa[e]:
                res_df[e].loc[c] = 2
                if res_df[c].loc[e] == 0:
                    res_df[c].loc[e] = 1
        return res_df


def granger_mv2(data, sig_level=0.05, maxlag=5, verbose=False):
    g = Granger(data, p=maxlag)
    res_df = g.fit(alpha=sig_level)
    return res_df


def generate_structure(N=1000):
    epsw = np.random.randn(N)**3
    epsx = np.random.randn(N)**3
    epsy = np.random.randn(N)**3

    x = np.zeros([N])
    y1 = np.zeros([N])
    w1 = np.zeros([N])
    y2 = np.zeros([N])
    w2 = np.zeros([N])

    for i in range(3,N):
        x[i] = 0.3*x[i-1]+0.5*epsx[i]
        y1[i] = 0.8*y1[i-1]+0.5*epsy[i]
        w1[i] = -0.6*w1[i-1]+0.8*y1[i-1]+0.8*x[i-2]+0.5*epsw[i]

        y2[i] = 0.8*y2[i-1]+0.5*epsy[i] + 0.8*x[i-1]
        w2[i] = -0.6*w2[i-1]+0.8*y2[i-1]+0.8*x[i-2]+0.5*epsw[i]

    x = pd.DataFrame(x, columns=["V1"])
    y1 = pd.DataFrame(y1, columns=["V2"])
    w1 = pd.DataFrame(w1, columns=["V3"])
    y2 = pd.DataFrame(y2, columns=["V2"])
    w2 = pd.DataFrame(w2, columns=["V3"])

    # data1 has a v structure and data2 has a chain structure
    data1 = pd.concat([x, y1, w1], axis=1, sort=False)
    data2 = pd.concat([x, y2, w2], axis=1, sort=False)
    return data1, data2




if __name__ == "__main__":
    print("simulated data")
    data1, data2 = generate_structure()
    # data1, _, _ = unfaithful_fork_linear_generator()


    print(data1)
    # G = Granger(data1)
    # res1 = G.fit(alpha=0.05)
    #
    # G = Granger(data2)
    # res2 = G.fit(alpha=0.05)

    res1 = granger_mv2(data1)

    res2 = granger_mv2(data2)


    print("Simulated data:")
    print(res1)
    print(res2)