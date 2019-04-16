import numpy as np
from scipy import stats


class univariateNormDistrib(object):

    def __init__(self, ):
        pass

    def uvn(self, x, miu, sigma):
        x_ = np.array(x).T
        # print(x_)
        pdf = 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(- (x_ - miu) ** 2 / (2 * sigma ** 2))
        # return np.reshape(pdf, np.array(pdf).shape[0])
        return pdf