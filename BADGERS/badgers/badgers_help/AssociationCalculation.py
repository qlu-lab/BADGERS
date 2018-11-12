import logging
import pandas
import numpy
import gzip
from numpy import dot as d
from timeit import default_timer as timer

from .. import Constants
from .. import Exceptions

from ..PredictionModel import  WDBQF

class ARF(object):
    """Association result format"""
    TRAIT = 0
    TRAIT_ID = 1
    ZSCORE = 2
    EFFECT_SIZE = 3
    SIGMA_G_2 = 4
    N_SNPS_IN_MODEL = 5
    N_SNPS_IN_COV = 6
    N_SNPS_USED = 7

    K_TRAIT = "trait"
    K_TRAIT_ID = "trait_id"
    K_ZSCORE = "zscore"
    K_EFFECT_SIZE = "effect_size"
    K_VAR_G = "var_g"
    K_N_SNPS_IN_MODEL = "n_snps_in_model"
    K_N_SNPS_IN_COV = "n_snps_in_cov"
    K_N_SNPS_USED = "n_snps_used"

    order=[(TRAIT,K_TRAIT),(TRAIT_ID,K_TRAIT_ID), (ZSCORE,K_ZSCORE), (EFFECT_SIZE,K_EFFECT_SIZE), (SIGMA_G_2, K_VAR_G), (N_SNPS_IN_MODEL, K_N_SNPS_IN_MODEL), (N_SNPS_IN_COV, K_N_SNPS_IN_COV), (N_SNPS_USED,K_N_SNPS_USED)]

class Context(object):
    def __init__(self): raise Exceptions.ReportableException("Tried to instantiate abstract context")
    def get_weights(self, trait): pass
    def get_covariance(self, trait, snps): pass
    def get_n_in_covariance(self, trait): pass
    def get_gwas(self, snps): pass
    def get_model_snps(self): pass
    def get_data_intersection(self): pass
    def provide_calculation(self, trait): pass
    def get_model_info(self): pass

def association(trait, context, return_snps=False):
    #capture context  
    n_snps_in_model, i, cov, snps,var,location,trait_id = context.provide_calculation(trait)
    snps_used = i[Constants.SNP]
    n_snps_used = len(snps_used)
    zscore, effect_size, sigma_g_2 = numpy.nan, numpy.nan, numpy.nan
    if n_snps_used > 0:
        i_weight = i[WDBQF.K_WEIGHT]
        i_zscore = i[Constants.ZSCORE]
        i_beta = i[Constants.BETA]
        # sigma from reference
        variances = cov
        snp_name = []
        x_value = []
        x_value_real = []
        i_sigma_l = numpy.sqrt(variances)
        fn = gzip.open(location ,"r")
        for lines in fn:
            tmp = lines.split()
            snp_name.append(tmp[0])
            value = []
            for kj in xrange(len(tmp)-1):
                value.append(float(tmp[kj+1]))
            x_value.append(value)
        for snps in snps_used:
            if snps in snp_name:
                index = snp_name.index(snps)
                x_value_real.append(x_value[index])
        array = numpy.array(x_value_real)
        ti = numpy.dot(array.transpose(),i_weight)
        var_2 = numpy.var(ti)

        sigma_g_2 = var_2
        if sigma_g_2 >0:
            try:
                zscore = numpy.sum(i_weight * i_zscore * i_sigma_l) / numpy.sqrt(sigma_g_2)
                effect_size = numpy.sum(i_weight * i_beta * (i_sigma_l**2))/ sigma_g_2

            except Exception as e:
                logging.log(9, "Unexpected exception when calculating zscore: %s, %s", trait, str(e))
    r = (trait,trait_id, zscore, effect_size, sigma_g_2, n_snps_in_model, 0, n_snps_used)  
    if return_snps:
        return r, set(snps_used)
    else:
        return r


def dataframe_from_results(results):
    if len(results) == 0:
        results = [[],[],[],[],[],[],[]]
    r = pandas.DataFrame({key: results[order] for order, key in ARF.order})
    r = r[[key for order,key in ARF.order]]
    return r
