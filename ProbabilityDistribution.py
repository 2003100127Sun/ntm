import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.stats import multivariate_normal
from sklearn.preprocessing import scale


class probabilityDistribution(object):
    
    def __init__(self, ):
        pass
    
    def multiNorm(self, miu):
        a = (pos - miu).dot(np.linalg.inv(sigma))  # np.linalg.inv()矩阵求逆
        b = np.expand_dims(pos - u, axis=3)
        Z = np.zeros((num, num), dtype=np.float32)
        for i in range(num):
            Z[i] = [np.dot(a[i, j], b[i, j]) for j in range(num)]
        Z = np.exp(Z * (-0.5)) / (2 * np.pi * (np.linalg.det(sigma)) ** (0.5))
        
    def univarNorm(self, x, miu, sigma, is_plot=False):
        x_ = np.array(x)
        # print(x_)
        pdf = np.exp(-(np.power(x_ - miu, 2)) / (2 * np.power(sigma, 2))) / (sigma * np.sqrt(2 * np.pi))
        if is_plot == True:
            bind = pd.DataFrame({'x': x_, 'y': np.array(pdf)})
            bind = bind.sort_values(by=['x'], ascending=True)
            # print(bind)
            plt.plot(bind['x'], bind['y'], color='grey', linestyle='solid')
            plt.hist(x_, bins=30, density='norm', facecolor='skyblue', ec="grey", alpha=0.6)
            plt.show()
        return pdf


if __name__ == "__main__":
    p = probabilityDistribution()
    x = scale([1.0333994708994844e-05, 7.10978835978841e-05, 5.643738977072275e-05, -1.5652557319224223e-05, 1.6837522045855344e-05, -7.220017636684433e-06, -6.60548941798942e-05, -2.9954805996472002e-05, 2.4994488536155403e-05, 3.30412257495596e-05, 7.0216049382716e-05, 3.508046737213475e-05, 4.0784832451499044e-05, 5.266203703703726e-05, -8.101851851851837e-06, 6.0736331569664665e-05, 4.681988536155195e-05, -1.4991181657848153e-05, -8.325066137566169e-05, 4.883156966490357e-05, 4.938271604938276e-05, -8.107363315696618e-05, -9.31437389770754e-06, 2.039241622574944e-05, 7.727072310405628e-05, 3.51906966490298e-05, -1.3695987654320616e-05, -1.4329805996473019e-06, 3.0065035273368825e-05, 4.103284832451554e-05, -1.8022486772486615e-05, 3.3399470899470333e-05, 1.1849647266309047e-06, -1.1078042328042389e-05, 1.581790123456918e-05, -7.798721340387265e-06, 3.8855820105820026e-05, 1.9565696649024553e-06, 2.0392416225748628e-05, 2.9431216931217005e-05, 3.7202380952381145e-05, 5.5362654320988114e-05, 9.138007054673735e-05, 5.4370590828925155e-05, -1.0141093474426953e-05, 7.864858906525519e-05, 1.4329805996472273e-05, -2.6372354497353602e-05, -4.6847442680776266e-05, -4.042658730158761e-05, 4.219025573192175e-05, -1.2400793650792143e-06, 1.2676366843034504e-05, 4.309964726631406e-05, 5.798059964726615e-05, 3.0836640211640064e-05, 4.819775132275137e-05, -1.5239197530863988e-05, 3.392305996472746e-05, -5.153218694885388e-05, 3.4171075837743144e-05, 4.106040564373815e-05, 5.660273368606676e-05, 7.090498236331539e-05, 6.429122574955892e-05, 3.5769400352733064e-05, -2.083333333333326e-05, 1.7388668430335405e-05, 6.148037918871275e-05, -2.700617283949797e-06, 6.980268959435579e-05, 2.0171957671956853e-05, 2.1852954144620782e-05, 3.130511463844769e-05, -8.258928571428512e-05, 5.577601410934716e-05, 4.783950617283963e-05, -3.0092592592592438e-05])
    mu = np.mean(x)
    sigma = np.std(x)
    # plt.hist(x, bins=100, density='norm')
    # plt.show()
    print(p.univarNorm(x, mu, sigma, is_plot=True))