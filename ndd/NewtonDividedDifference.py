import numpy as np


class newtonDividedDifference(object):

    def __init__(self, ):
        pass

    def getDivDiffCrafted(self, x, y):
        """
        get lower Triangular Matrix for Divided Difference Table.
        | -
        |     -
        |         -
        |             -
        |----------------
        .. @Source:
           --------
           http://dmpeli.math.mcmaster.ca/Matlab/Math1J03/LectureNotes/Lecture4_1.htm
           & answer 1 from
           https://stackoverflow.com/questions/14823891/newton-s-interpolating-polynomial-python.

        :param x: 1d list
        :param y: 1d list
        :return: 2d array
        """
        x_ = np.array(x)
        y_ = np.array(y)
        n = x_.shape[0]
        D = np.zeros([n, n])
        D[:, 0] = np.transpose(y_)
        for i in range(1, n):
            for j in range(i, n):
                # print(j)
                D[j, i] = (D[j, i-1] - D[j-1, i-1]) / (x_[j] - x_[j-i])
                # print(D[j, i])
        return D
        
    def getDivDiffTable(self, x, y):
        """
        get upper Triangular Matrix for Standard Divided Difference Table in accordance
        with the table form in 'Numerical Analysis' by Zhizhong Sun.

        .. @Example:
           ---------
           >>> x = [5, 6, 9, 11]
           >>> y = [12, 13, 14, 16]
           >>> D[1, 4] = 13.47

        .. @Source:
           --------
           https://www.geeksforgeeks.org/newtons-divided-difference-interpolation-formula/

        :param x: 1d list
        :param y: 1d list
        :return: 2d array (Standard Divided Difference Table)
        """
        x_ = np.array(x)
        y_ = np.array(y)
        n = x_.shape[0]
        D = np.zeros([n, n])
        D[:, 0] = np.transpose(y_)
        for i in range(1, n):
            for j in range(n-i):
                D[j][i] = ((D[j][i-1] - D[j+1][i-1]) / (x_[j] - x_[i+j]))
        return D

    def getDivDiffDiag(self, x, y):
        """
        get diagonal elements from a divided difference table.
        ----------------
        |            -
        |         -
        |      -
        |   -
        |-

        .. @Source:
           --------
           https://stackoverflow.com/questions/14823891/newton-s-interpolating-polynomial-python

        :param x: 1d list
        :param y: 1d list
        :return: list.
        """
        x_ = np.copy(np.array(x))
        n = x_.shape[0]
        y_ = np.copy(np.array(y))
        for k in range(1, n):
            y_[k:n] = (y_[k:n] - y_[k-1]) / (x_[k:n] - x_[k-1])
        return y_


if '__main__' in __name__:
    p = newtonDividedDifference()
    x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    # y = [1, 1.5, 2, 6, 10.5, 11, 15, 31, 33, 34]
    # y = [0.091122961, 0.06913842, 0.052153563, 0.039165723, 0.029312231, 0.021881271, 0.016302499, 0.012128435, 0.009013299, 0.006692851]
    y = [4, 7.1, 10.8, 15, 21.6, 25, 29, 31, 33, 34]
    print(p.getDivDiffCrafted(x, y))
    # print(p.getDivDiffDiag(x, y))
    # print(p.getDivDiffTable(x, y))