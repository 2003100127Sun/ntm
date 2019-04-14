import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import rpy2.robjects as rob
from ntim.bcNTMapping import bcNTMapping


class batchPermutation(object):

    def __init__(self, ):
        self.bc_nt = bcNTMapping()

    def qtl(self, ):
        return rob.r(
            '''
            poplar <- read.cross(
                "csvs",
                ".",
                geno,
                pheno,
                genotypes=c("AA", "AB")
            )
            gp <- calc.genoprob(poplar, step=1, error.prob=0.001)
            p.em <- scanone(gp)
            p.hk <- scanone(gp, method="hk")
            p.ehk <- scanone(gp, method="ehk")
            # plot(p.em, p.hk, p.ehk, ylab="LOD score", lty=c(1,1,2))
            # plot(p.em, p.hk, ylab="LOD score", lty=c(1,1,2))
            return ( list('em' = p.em, 'hk' = p.hk, 'ehk' = p.ehk) )
            '''
        )

    def plot_(self, i, flag, em=None, hk=None):
        if em is not None:
            plt.plot(em['pos'], em['lod'], c='coral', lw=3, alpha=0.75, label='Standard IM')
            plt.xlabel('Genomic distance (cM)', fontsize=12)
            plt.ylabel('LOD', fontsize=12)
            if flag == 1:
                plt.title('Group' + str(12 + 1) + ' with original data')
            if flag == 2:
                plt.title('Group' + str(12 + 1) + ' with NT measurement')
        if hk is not None:
            plt.plot(hk['pos'], hk['lod'], c='b', label='Haley-Knott regression')
        plt.legend(fontsize=10)
        plt.show()

    def probability(self, num_groups_chr, num_permutation=500):
        # for i in range(num_permutation):
        #     self.bc_nt.enhanced(sv_fpn='./data/populus.BC.pheno_e' + str(i+1) + '.csv')
        rob.r(
            '''
            library(qtl)
            setwd('./data/')
            '''
        )
        accs = []
        em_lods = []
        for i in range(num_groups_chr):
            print('Group {}:'.format(i+1))
            dist = bcNTMapping(
                pheno_data='populus.BC.pheno_ori78.csv',
                geno_data="populus.BC.geno" + str(i+1) + ".csv"
            ).genomic_dists
            print('genomic dist is: {}'.format(dist))
            rob.globalenv["geno"] = "populus.BC.geno" + str(i+1) + ".csv"
            rob.globalenv["pheno"] = "populus.BC.pheno_ori78.csv"
            qtl_analysis = self.qtl()
            em_marks = rob.r.rownames(qtl_analysis[0])
            em_group = qtl_analysis[0][0]
            em_pos = qtl_analysis[0][1]
            em_lod = qtl_analysis[0][2]
            # print(em_marks)
            hk_marks = rob.r.rownames(qtl_analysis[1])
            hk_group = qtl_analysis[1][0]
            hk_pos = qtl_analysis[1][1]
            hk_lod = qtl_analysis[1][2]
            em = pd.DataFrame({
                'group': np.array(em_group),
                'pos': np.array(em_pos),
                'lod': np.array(em_lod)
            }, index=em_marks
            )
            hk = pd.DataFrame({
                'group': np.array(hk_group),
                'pos': np.array(hk_pos),
                'lod': np.array(hk_lod)
            }, index=hk_marks
            )
            self.plot_(i, 1, em=em)
            em_id_max = em['lod'].idxmax()
            em_pos_max = em.loc[em_id_max, 'pos']
            em_lod_max = em['lod'].max()
            # print(em_nt_max)
            # print(em_id)
            conf_interv = []
            print('em_pos_max: {}'.format(em_pos_max))
            for k in range(len(dist)):
                if em_pos_max == dist[k]:
                    print(dist[k])
                    if k == 0 and k != len(dist) - 1:
                        conf_interv = [dist[k], dist[k + 1]]
                    elif k != 0 and k == len(dist) - 1:
                        conf_interv = [dist[k-1], dist[k]]
                    else:
                        conf_interv = [dist[k-1], dist[k+1]]
                else:
                    if em_pos_max > dist[k] and em_pos_max < dist[k+1]:
                        conf_interv = [dist[k], dist[k+1]]
            print('confident interval: {}'.format(conf_interv))
            arg_max = []
            for j in range(num_permutation):
                rob.globalenv["pheno"] = "populus.BC.pheno_e" + str(j+1) + '.csv'
                qtl_analysis = self.qtl()
                em_marks = rob.r.rownames(qtl_analysis[0])
                em_group = qtl_analysis[0][0]
                em_pos = qtl_analysis[0][1]
                em_lod = qtl_analysis[0][2]
                hk_marks = rob.r.rownames(qtl_analysis[1])
                hk_group = qtl_analysis[1][0]
                hk_pos = qtl_analysis[1][1]
                hk_lod = qtl_analysis[1][2]
                em = pd.DataFrame({
                    'group': np.array(em_group),
                    'pos': np.array(em_pos),
                    'lod': np.array(em_lod)
                }, index=em_marks
                )
                hk = pd.DataFrame({
                    'group': np.array(hk_group),
                    'pos': np.array(hk_pos),
                    'lod': np.array(hk_lod)
                }, index=hk_marks
                )
                em_nt_id_max = em['lod'].idxmax()
                em_nt_pos_max = em.loc[em_nt_id_max, 'pos']
                print('em_nt_pos_max: {}'.format(em_nt_pos_max))
                # em_nt_max = em['lod'].max()
                if conf_interv[0] <= em_nt_pos_max <= conf_interv[1]:
                    arg_max.append(1)
                else:
                    arg_max.append(0)
                self.plot_(i, 2, em=em)
                # print(em_marks)
            acc = sum(arg_max) / len(arg_max)
            print(acc)
            em_lods.append([em_id_max, em_pos_max, em_lod_max, acc])
        # pd.DataFrame(em_lods).to_csv('ex4.nt', sep='\t', header=None, index=False)
        return accs


if __name__ == "__main__":
    p = batchPermutation()
    print(p.probability(num_groups_chr=22, num_permutation=1000))