import numpy as np
import pandas as pd


class mapping(object):

    def __init__(self, ):
        pass

    def read(self, file_path, sep=','):
        return pd.read_csv(file_path, sep=sep, header=0)

    def phenotype(self, pheno_data, sep=','):
        phenotypic_data = self.read(pheno_data, sep)
        drop_na = phenotypic_data.dropna(how='any')
        indices_na = phenotypic_data.loc[phenotypic_data.isnull().any(axis=1) == True].index
        drop_na = drop_na.set_index(np.arange(drop_na.shape[0]), inplace=False, drop=True)
        return drop_na, indices_na

    def genotype(self, geno_data, indices_na, sep=','):
        """

        :param geno_data:
        :param indices_na:
        :param sep:
        :return:
        """
        genotypic_data = self.read(geno_data, sep)
        chromosome_groups = genotypic_data.iloc[0, 1:]
        genomic_dists = genotypic_data.iloc[1, 1:]
        if len(indices_na):
            drop_na = genotypic_data.drop(index=indices_na + 2)
            drop_na = drop_na.drop(index=[0, 1])
        else:
            drop_na = genotypic_data.drop(index=[0, 1])
        drop_na = drop_na.set_index(np.arange(drop_na.shape[0]), inplace=False, drop=True)
        return drop_na, chromosome_groups.astype(np.int64).tolist(), genomic_dists.astype(np.float64).tolist()