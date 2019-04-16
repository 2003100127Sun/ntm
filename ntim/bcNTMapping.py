import sys
sys.path.append('../')
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ndd.NewtonDividedDifference import newtonDividedDifference
from distributions.UnivariateNormDistrib import univariateNormDistrib
from sklearn.preprocessing import scale
from sklearn import mixture
from ntim.base.Mapping import mapping


class bcNTMapping(mapping):

    def __init__(self, pheno_data='./populus.BC.pheno.csv', geno_data='./populus.BC.geno.csv'):
        super(bcNTMapping, self).__init__()
        self.parser = argparse.ArgumentParser(description='Sequence Identity calculations')
        self.parser.add_argument(
            "-fp", "--filepheno", metavar='file_pheno',
            dest='fp', help='filepheno', type=str
        )
        self.parser.add_argument(
            "-fg", "--filegeno", metavar='file_geno',
            dest='fg', help='filegeno', type=str
        )
        self.parser.add_argument(
            "-nt", "--ntvalues", metavar='nt_values',
            dest='nt', help='numerical trajectory values', type=int
        )
        self.parser.add_argument(
            "-spt", "--split", metavar='split',
            dest='spt', help='split', type=int
        )
        self.parser.add_argument(
            "-it", "--iTables", metavar='interval_table',
            dest='it', help='interval table', type=int
        )
        args = self.parser.parse_args()
        if args.fp:
            self.pheno_data = args.fp
        else:
            self.pheno_data = pheno_data
        if args.fg:
            self.geno_data = args.fg
        else:
            self.geno_data = geno_data
        self.und = univariateNormDistrib()
        self.phenos, self.indices_na = self.phenotype(pheno_data=self.pheno_data)
        self.genos, self.chromosome_groups, self.genomic_dists = self.genotype(geno_data=self.geno_data, indices_na=self.indices_na)
        self.num_samples = self.phenos.shape[0]
        self.num_intervs = self.genos.shape[1] - 2
        self.num_genos_qtl = 2
        if args.nt == 1:
            print('NT values using NT measurement: \n{}'.format(self.ntvs()))
        if args.spt == 1:
            print('Spliting...: {}'.format(self.split()))
        if args.it == 1:
            print('Intervals: {}'.format(self.intervals()))
            print('Interval table: {}'.format(self.iTable()))

    def split(self, ):
        previous = 0
        num = self.num_intervs + 1
        genos = self.genos.iloc[:, 1:]
        id = pd.DataFrame(np.arange(genos.shape[0]) + 1, columns=['id'])
        id.loc[-2] = ''
        id.loc[-1] = ''
        for i in range(num):
            if i+1 < num and self.genomic_dists[i+1] - self.genomic_dists[i] < 0:
                chr = self.chromosome_groups[previous: i + 1]
                dist = self.genomic_dists[previous: i + 1]
                genos_chr = genos.iloc[:, previous: i + 1]
                genos_chr = genos_chr.astype(np.str)
                # print(chr)
                # print(dist)
                for j in range(genos_chr.shape[0]):
                    for k in range(genos_chr.shape[1]):
                        if genos_chr.iloc[j, k] == '-1':
                            genos_chr.loc[j][k] = '-'
                        elif genos_chr.iloc[j, k] == '1':
                            genos_chr.loc[j][k] = 'AB'
                        elif genos_chr.iloc[j, k] == '0':
                            genos_chr.loc[j][k] = 'AA'

                genos_chr.loc[-2] = chr
                genos_chr.loc[-1] = dist
                genos_chr = genos_chr.sort_index()
                s = pd.concat([id, genos_chr], axis=1)
                pd.DataFrame(s).to_csv('./populus.BC.geno' + str(chr[0]) + '.csv', sep=',', header=True, index=None)
                previous = i + 1
        chr = self.chromosome_groups[previous: num + 1]
        dist = self.genomic_dists[previous: num + 1]
        last_chr = genos.iloc[:, previous: num + 1]
        last_chr = last_chr.astype(np.str)
        for j in range(last_chr.shape[0]):
            for k in range(last_chr.shape[1]):
                if last_chr.iloc[j, k] == '-1':
                    last_chr.loc[j][k] = '-'
                elif last_chr.iloc[j, k] == '1':
                    last_chr.loc[j][k] = 'AB'
                elif last_chr.iloc[j, k] == '0':
                    last_chr.loc[j][k] = 'AA'
        print(last_chr)
        last_chr.loc[-2] = chr
        last_chr.loc[-1] = dist
        last_chr = last_chr.sort_index()
        s = pd.concat([id, last_chr], axis=1)
        pd.DataFrame(s).to_csv('./populus.BC.geno' + str(chr[0]) + '.csv', sep=',', header=True, index=None)
        return 0

    def ntvs(self):
        x = np.arange(1, self.phenos.shape[1])
        div_diff = []
        for i in range(self.phenos.shape[0]):
            dd = newtonDividedDifference().getDivDiffDiag(x, np.array(self.phenos.iloc[i, 1:]))
            div_diff.append(dd[-1])
            # print(dd[-1])
        div_diff = np.array(div_diff)[:, np.newaxis]
        return div_diff

    def intervals(self):
        marks = self.genos.iloc[:, 1:]
        intervals = {}
        for i in range(self.num_intervs):
            intervals[i+1] = np.array(marks.iloc[:, i:i+2])
        return intervals

    def iTable(self):
        intervs = self.intervals()
        interv_Table = np.zeros([self.num_samples, self.num_intervs]) - 1
        for i in range(self.num_intervs):
            for j in range(self.num_samples):
                if intervs[i+1][j][0] == 1 and intervs[i+1][j][1] == 1:
                    interv_Table[j][i] = 1
                elif intervs[i+1][j][0] == 1 and intervs[i+1][j][1] == 0:
                    interv_Table[j][i] = 2
                elif intervs[i+1][j][0] == 0 and intervs[i+1][j][1] == 1:
                    interv_Table[j][i] = 3
                elif intervs[i+1][j][0] == 0 and intervs[i+1][j][1] == 0:
                    interv_Table[j][i] = 4
                else:
                    interv_Table[j][i] = 5
        return interv_Table.astype(np.int64)
    
    def yTable(self, interv_cls, phenos_enhanced):
        y_table = {}
        phenos_ = np.squeeze(phenos_enhanced, axis=(1, ))
        pheno1 = []
        pheno2 = []
        pheno3 = []
        pheno4 = []
        for i in range(self.num_samples):
            if interv_cls[i] == 1:
                pheno1.append(phenos_[i])
            elif interv_cls[i] == 2:
                pheno2.append(phenos_[i])
            elif interv_cls[i] == 3:
                pheno3.append(phenos_[i])
            elif interv_cls[i] == 4:
                pheno4.append(phenos_[i])
        y_table['pheno1'] = pheno1
        y_table['pheno2'] = pheno2
        y_table['pheno3'] = pheno3
        y_table['pheno4'] = pheno4
        return y_table

    def llreduced(self, phenos_enhanced, miu, sigma):
        phenos_ = np.squeeze(phenos_enhanced, axis=(1,))
        return np.sum(np.log(self.und.uvn(x=phenos_, miu=miu, sigma=sigma)))
    
    def Pij1(self, phenos, miu1, miu2, sigma, theta):
        prop1 = (1 - theta) * self.und.uvn(x=phenos, miu=miu1, sigma=sigma)
        prop2 = theta * self.und.uvn(x=phenos, miu=miu2, sigma=sigma)
        return prop1 / (prop1 + prop2)

    def Pij2(self, phenos, miu1, miu2, sigma, theta):
        prop1 = theta * self.und.uvn(x=phenos, miu=miu1, sigma=sigma)
        prop2 = (1 - theta) * self.und.uvn(x=phenos, miu=miu2, sigma=sigma)
        return prop1 / (prop1 + prop2)

    def mstep(self, pheno2, pheno3, miu1, miu2, sigma, theta):
        gamma = [0 for _ in range(4)]
        gamma[0] = self.Pij1(pheno2, miu1, miu2, sigma, theta)
        gamma[1] = 1 - gamma[0]
        gamma[2] = self.Pij2(pheno3, miu1, miu2, sigma, theta)
        gamma[3] = 1 - gamma[2]
        return gamma

    def estep(self, pheno1, pheno2, pheno3, pheno4, gamma):
        if len(pheno1):
            n1 = len(pheno1)
        else:
            n1 = 0
        if len(pheno2):
            n2 = len(pheno2)
        else:
            n2 = 0
        if len(pheno3):
            n3 = len(pheno3)
        else:
            n3 = 0
        if len(pheno4):
            n4 = len(pheno4)
        else:
            n4 = 0
        n = n1 + n2 + n3 + n4
        miu1_ = (np.sum(pheno1) + np.sum(gamma[0] * pheno2) + np.sum(gamma[2] * pheno3)) / (
                    n1 + np.sum(gamma[0]) + np.sum(gamma[2]))
        miu2_ = (np.sum((gamma[1]) * pheno2) + np.sum((gamma[3]) * pheno3) + np.sum(pheno4)) / (
                np.sum((gamma[1])) + np.sum((gamma[3])) + n4)
        part1 = np.sum(np.power((pheno1 - miu1_), 2))
        part2 = np.sum(gamma[0] * np.power((pheno2 - miu1_), 2) + (gamma[1]) * np.power((pheno2 - miu2_), 2))
        part3 = np.sum(gamma[2] * np.power((pheno3 - miu1_), 2) + (gamma[3]) * np.power((pheno3 - miu2_), 2))
        part4 = np.sum(np.power((pheno4 - miu2_), 2))
        sigma_ = np.sqrt((part1 + part2 + part3 + part4) / n)
        theta_ = (np.sum(gamma[1]) + np.sum(gamma[2])) / (n2 + n3)
        return miu1_, miu2_, sigma_, theta_

    def w(self, theta):
        w = dict()
        w[0][0] = 1
        w[0][1] = 0
        w[1][0] = 1 - theta
        w[1][1] = theta
        w[2][0] = theta
        w[2][1] = 1 - theta
        w[3][0] = 0
        w[3][1] = 1
        # for i in range(self.num_samples):
        #     for j in range(self.num_genos_qtl):
        #         w[i][j] =
        return w

    def mle(self, pheno1, pheno2, pheno3, pheno4, miu1, miu2, sigma, theta):
        """
        mle consists of m step and e step.
        .. @past
           -----
            # print('pheno lens: {}, {}, {}, {}'.format(n1, n2, n3, n4))
            # w1 = self.Pij1(pheno2, miu1, miu2, sigma, theta)
            # w2 = self.Pij2(pheno3, miu1, miu2, sigma, theta)
            # miu1_ = (np.sum(pheno1) + np.sum(w1 * pheno2) + np.sum(w2 * pheno3)) / (n1 + np.sum(w1) + np.sum(w2))
            # miu2_ = (np.sum((1 - w1) * pheno2) + np.sum((1 - w2) * pheno3) + np.sum(pheno4)) / (np.sum((1 - w1)) + np.sum((1 - w2)) + n4)
            # part1 = np.sum(np.power((pheno1 - miu1_), 2))
            # part2 = np.sum(w1 * np.power((pheno2 - miu1_), 2) + (1 - w1) * np.power((pheno2 - miu2_), 2))
            # part3 = np.sum(w2 * np.power((pheno3 - miu1_), 2) + (1 - w2) * np.power((pheno3 - miu2_), 2))
            # part4 = np.sum(np.power((pheno4 - miu2_), 2))
            # sigma_ = np.sqrt((part1 + part2 + part3 + part4) / n)
            # theta_ = (np.sum(1 - w1) + np.sum(w2)) / (n2 + n3)

        :param pheno1:
        :param pheno2:
        :param pheno3:
        :param pheno4:
        :param miu1:
        :param miu2:
        :param sigma:
        :param theta:
        :return:
        """
        ### /* M Step */ ###
        gamma = self.mstep(pheno2, pheno3, miu1, miu2, sigma, theta)
        ### /* E Step */ ###
        miu1_, miu2_, sigma_, theta_ = self.estep(pheno1, pheno2, pheno3, pheno4, gamma)
        return miu1_, miu2_, sigma_, theta_

    def logL(self, pheno1, pheno2, pheno3, pheno4, miu1, miu2, sigma, theta):
        part1 = np.sum(np.log(self.und.uvn(x=pheno1, miu=miu1, sigma=sigma)))
        part2 = np.sum(
            (1 - theta) * np.log(self.und.uvn(x=pheno2, miu=miu1, sigma=sigma)) + theta * np.log(self.und.uvn(x=pheno2, miu=miu2, sigma=sigma))
        )
        part3 = np.sum(
            theta * np.log(self.und.uvn(x=pheno3, miu=miu1, sigma=sigma)) + (1 - theta) * np.log(self.und.uvn(x=pheno3, miu=miu2, sigma=sigma))
        )
        part4 = np.sum(np.log(self.und.uvn(x=pheno4, miu=miu2, sigma=sigma)))
        # print('logL part1: {}'.format(part1))
        # print('logL part2: {}'.format(part2))
        # print('logL part3: {}'.format(part3))
        # print('logL part4: {}'.format(part4))
        return part1 + part2 + part3 + part4

    def lod(self, reduced, full):
        """
        :param reduced: numerator of mle
        :param full: denominator of mle
        :return: float
        """
        return 2 * (full - reduced)
    
    def localize(self, interv, theta):
        dist_intervs = {}
        for i in range(self.num_intervs):
            dist_intervs[i + 1] = self.genomic_dists[i:i + 2]
        rf = 0.5 * (1 - np.exp(-2 * (dist_intervs[interv+1] - dist_intervs[interv])))
        localized_dist = -0.5 * np.log(1 - 2 * rf * theta)
        return localized_dist

    def predict(self, num_epochs=250, num_optimize=100, ntim=True, rand=True):
        itable = p.iTable()
        LRT = []
        LRs = []
        for i in range(num_epochs):
            LR = []
            if ntim and rand:
                pheno_sampling = scale(np.array(pd.DataFrame(self.ntvs()).sample(self.num_samples)))
            elif ntim is True and rand is False:
                pheno_sampling = np.abs(scale(np.array(self.ntvs())))
            elif ntim is False and rand is False:
                pheno_sampling = np.array(pd.DataFrame(self.phenos.iloc[:, 11]))
            else:
                pheno_sampling = np.array(pd.DataFrame(self.phenos.iloc[:, -1]).sample(self.num_samples))
            for j in range(self.num_intervs):
                pheno1 = p.yTable(itable[:, j], pheno_sampling)['pheno1']
                pheno2 = p.yTable(itable[:, j], pheno_sampling)['pheno2']
                pheno3 = p.yTable(itable[:, j], pheno_sampling)['pheno3']
                pheno4 = p.yTable(itable[:, j], pheno_sampling)['pheno4']
                # print(pheno1, pheno2, pheno3, pheno4)
                # if len(pheno1):
                #     pheno1 = np.array(pheno1)
                # else:
                #     pheno1 = np.nan
                # if len(pheno2):
                #     pheno2 = pheno2
                # else:
                #     pheno2 = np.nan
                # if len(pheno3):
                #     pheno3 = pheno3
                # else:
                #     pheno3 = np.nan
                # if len(pheno4):
                #     pheno4 = pheno4
                # else:
                #     pheno4 = np.nan
                miu1 = np.mean(pheno1)
                miu2 = np.mean(pheno4)
                miu = np.mean(pheno_sampling)
                sd = np.std(pheno_sampling)
                theta = 0.25
                for k in range(num_optimize):
                    miu1, miu2, sd, theta = self.mle(np.array(pheno1), np.array(pheno2), np.array(pheno3), np.array(pheno4), miu1, miu2, sd, theta)
                    # print(miu1, miu2, sd, theta)
                LR.append(
                    self.lod(
                        full=self.logL(pheno1, pheno2, pheno3, pheno4, miu1, miu2, sd, theta),
                        reduced=self.llreduced(pheno_sampling, miu, sd)
                    )
                )
            print(LR)
            LRs.append(LR)
            LRT.append(np.max(LR))
        return LRs, LRT

    def plots(self, ):
        # dd = [13.711407689738195, 22.55570603331745, 30.556589086018306, 28.77699416074026, 12.725182435636668, 24.588766699055952, 24.57098811370858]
        dd = [5.784548638030259, 10.795167866055209, 21.54282783231706, 20.140576866963514, 8.525692149076804, 16.544543239569435, 13.057609897626293]
        # plt.plot(np.arange(8), dd[118:126])
        plt.plot(np.arange(7), dd)
        print(dd)
        plt.show()

    def labels(self, ):
        gmm = mixture.GaussianMixture(n_components=2, covariance_type='full')
        gmm.fit(scale(self.ntvs()))
        pred_lables = gmm.predict(scale(self.ntvs()))
        ori = gmm.fit(self.phenos.iloc[:, -1][:, np.newaxis])
        ori_miu = ori.means_
        ori_sigma = ori.covariances_
        ori_lables = gmm.predict(self.phenos.iloc[:, -1][:, np.newaxis])
        # ori_lab = gmm.predict_proba(self.phenos.iloc[:, -1][:, np.newaxis])
        # print(ori_lables)
        # print(pred_lables)
        # print(pred0, pred1)
        ori0 = ori_lables[ori_lables == 0].shape[0]
        ori1 = ori_lables[ori_lables == 1].shape[0]
        pred0 = pred_lables[pred_lables == 0].shape[0]
        pred1 = pred_lables[pred_lables == 1].shape[0]
        if ori1 > ori0 and pred1 < pred0:
            pred_lables = np.array([0 if pred_lables.astype(np.bool)[i] else 1 for i in range(pred_lables.shape[0])])
        elif ori1 < ori0 and pred1 > pred0:
            pred_lables = np.array([0 if pred_lables.astype(np.bool)[i] else 1 for i in range(pred_lables.shape[0])])
        # print(pred_lables)
        # print(ori0, ori1)
        # print(pred0, pred1)
        return pred_lables, ori_lables, ori_miu, ori_sigma

    def enhanced(self, sv_fpn):
        pred_l, ori_l, ori_miu, ori_sigma = self.labels()
        pheno_mature = self.phenos.iloc[:, -1][:, np.newaxis]
        pheno_enhanced = []
        beta = np.random.randn()
        for i in range(self.phenos.shape[0]):
            if ori_l[i] == 0:
                if ori_l[i] == pred_l[i]:
                    pheno_enhanced.append(pheno_mature[i][0])
                else:
                    # beta = np.random.randn()
                    alpha = beta * (ori_miu[1] - ori_miu[0]) / (pheno_mature[i][0] - ori_miu[0]) + 1
                    enhanced = ori_miu[0] + alpha * (pheno_mature[i][0] - ori_miu[0])
                    pheno_enhanced.append(np.around(enhanced.item(0), 1))
            else:
                if ori_l[i] == pred_l[i]:
                    pheno_enhanced.append(pheno_mature[i][0])
                else:
                    # beta = np.random.randn()
                    alpha = beta * (ori_miu[0] - ori_miu[1]) / (pheno_mature[i][0] - ori_miu[1]) + 1
                    enhanced = ori_miu[1] + alpha * (pheno_mature[i][0] - ori_miu[1])
                    pheno_enhanced.append(np.around(enhanced.item(0), 1))
        # print(ori_miu)
        # print(pheno_enhanced)
        id = np.arange(len(pheno_enhanced)) + 1
        pd.DataFrame({'pheno': pheno_enhanced, 'id': id}).to_csv(sv_fpn, sep=',', header=True, index=False)
        return


if __name__ == "__main__":
    p = bcNTMapping(
        # pheno_data='./populus.BC.pheno.csv',
        # geno_data='./populus.BC.geno.csv'

        pheno_data='./poplar.pheno.csv',
        geno_data='./data/populus.BC.geno1.csv'
    )
    # 'python bcNTMapping.py -fp ./poplar.pheno.csv -fg ./data/populus.BC.geno1.csv'

    # print(p.genos)

    # print(p.chromosome_groups)

    print('genomic distance: {}'.format(p.genomic_dists))

    # print(p.ntvs())

    # print(p.split())

    # print(p.intervals())

    # print(p.phenos)

    # print(p.phenos.iloc[:, 11])

    # print(pd.DataFrame(np.abs(p.ntvs()*1000000)).to_csv('ss.csv', sep=',', header=None))

    # print(p.iTable())

    # print(p.localize())

    # itable = p.iTable()
    # print('itable:\n {}'.format(itable))
    # pheno1 = p.yTable(itable[:, 0], p.ntvs())['pheno1']
    # pheno2 = p.yTable(itable[:, 0], p.ntvs())['pheno2']
    # pheno3 = p.yTable(itable[:, 0], p.ntvs())['pheno3']
    # pheno4 = p.yTable(itable[:, 0], p.ntvs())['pheno4']
    # print('pheno lens: {}, {}, {}, {}'.format(len(pheno1), len(pheno2), len(pheno3), len(pheno4)))
    # print('phenos partition: {}, {}, {}, {}'.format(pheno1, pheno2, pheno3, pheno4))
    # # print(itable[:, 0][itable[:, 0] == 4])
    # print(p.mle(pheno1, pheno2, pheno3, pheno4, 19, 19.3, 4, 0.3))
    # print(p.logL(pheno1, pheno2, pheno3, pheno4, 19, 19.3, 4, 0.3))

    # print(p.predict(num_epochs=250, num_optimize=100, ntim=True, rand=False))

    # print(p.plots())

    # l1, l2, _, _ = p.labels()
    # print(l1[l1 == 1].shape)
    # print(l2[l2 == 1].shape)

    # print(p.enhanced(
    #     sv_fpn='./r/populus.BC.pheno_e.csv'
    # ))

    # print(p.probability(num_permutation=20))