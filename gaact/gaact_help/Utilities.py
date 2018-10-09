#author:'Donghui'

import logging
import pandas
import os
import gzip
import numpy
from scipy import stats

from .. import Constants
from .. import Utilities
from .. import MatrixManager
from ..PredictionModel import WDBQF, load_model, dataframe_from_weight_data

import AssociationCalculation


class OptimizedContext:
    def __init__(self, model_db_path,cov_loc,gwas, model,covariance):
        self.model_db_path = model_db_path
        self.cov_loc = cov_loc
        self.weight_data, self.snps_in_model = _prepare_weight_data(model)
        self.gwas_data = _prepare_gwas_data(gwas)
        self.covariance = covariance

    def get_covariance(self, trait, snps):
        return self.covariance.get(trait, snps, strict=False)

    def _get_weights(self, trait):
        w = self.weight_data[trait]
        w = {x[WDBQF.RSID]:x[WDBQF.WEIGHT] for x in w}
        return w

    def get_weights(self, trait):
        w = self.weight_data[trait]
        w = dataframe_from_weight_data(zip(*w))
        return w

    def get_model_snps(self):
        return set(self.snps_in_model)

    def _get_gwas(self, snps):
        snps = set(snps)
        g = self.gwas_data
        g = [g[x] for x in snps if x in g]
        g = {x[0]:(x[1], x[2]) for x in g}
        return g

    def get_gwas(self, snps):
        snps = set(snps)
        g = self.gwas_data
        g = [g[x] for x in snps if x in g]
        if len(g):
            g = zip(*g)
            g = pandas.DataFrame({Constants.SNP:g[0], Constants.ZSCORE:g[1], Constants.BETA:g[2]})
        else:
            g = pandas.DataFrame(columns=[Constants.SNP, Constants.ZSCORE, Constants.BETA])
        return g

    def get_data_intersection(self):
        return _data_intersection_2(self.weight_data, self.gwas_data)

    def provide_calculation(self, trait):
        w = self._get_weights('all')
        gwas = self._get_gwas(w.keys())
        type = [numpy.str, numpy.float64, numpy.float64, numpy.float64]
        columns = [Constants.SNP, WDBQF.K_WEIGHT, Constants.ZSCORE, Constants.BETA]
        d = {x: v for x, v in w.iteritems() if x in gwas}
        snps, cov = self.get_covariance(trait, d.keys())
        if snps is None:
            d = pandas.DataFrame(columns=columns)
            return len(w), d, cov, snps

        d = [(x, w[x], gwas[x][0], gwas[x][1]) for x in snps]
        d = zip(*d)
        head, tail = os.path.split(self.cov_loc)
        file_name = tail.split(".")[0]
        file_id = ("_").join(file_name.split("_")[:-1])
        kn = head+"/"+file_name+"_tn.txt.gz"
        #fw = open(head+"/"+file_name+"_tk.txt","r")
        #for line in fw:
        #        tmp = line.split(",")
        #        for kk in tmp:
        #            var.append(float(kk))
        #fw.close();
        #var = numpy.var(var)
        #print(var)
        var = 0
        if len(d):
            d = {columns[i]:numpy.array(d[i], dtype=type[i]) for i in xrange(0,len(columns))}
        else:
            d = {columns[i]:numpy.array([]) for i in xrange(0,len(columns))}
        return  len(w), d, cov, snps,var,kn,file_id



def _data_intersection_2(weight_data, gwas_data):
    traits = set()
    snps = set()
    for trait, entries in weight_data.iteritems():
        gs = zip(*entries)[WDBQF.RSID]
        for s in gs:
            if s in gwas_data:
                traits.add(trait)
                snps.add(s)
    return traits, snps

def _sanitized_gwas(gwas):
    gwas = gwas[[Constants.SNP, Constants.ZSCORE, Constants.BETA]]
    if numpy.any(~ numpy.isfinite(gwas[Constants.ZSCORE])):
        logging.warning("Discarding non finite GWAS zscores")
        gwas = gwas.loc[numpy.isfinite(gwas[Constants.ZSCORE])]
    return gwas

def _prepare_gwas(gwas):
    #If zscore is numeric, then everything is fine with us.
    # if not, try to remove "NA" strings.
    try:
        i = gwas.zscore.apply(lambda x: x != "NA")
        gwas = gwas.loc[i]
        gwas = pandas.DataFrame(gwas)
        gwas.loc[:,Constants.ZSCORE] = gwas.zscore.astype(numpy.float64)
    except Exception as e:
        logging.info("Unexpected issue preparing gwas... %s", str(e))
        pass

    if not Constants.BETA in gwas:
        gwas.loc[:,Constants.BETA] = numpy.nan

    return gwas

def _prepare_gwas_data(gwas):
    data = {}
    for x in gwas.values:
        data[x[0]] = x
    return data

def _prepare_model(model):
    K = WDBQF.K_TRAIT
    g = model.weights[K]
    model.weights[K] = pandas.Categorical(g, g.drop_duplicates())
    return model

def _prepare_weight_data(model):
    d = {}
    snps = set()
    for x in model.weights.values:
        trait = "all"
        if not trait in d:
            d[trait] = []
        entries = d[trait]
        entries.append(x)
        snps.add(x[WDBQF.RSID])
    return d, snps

def _beta_loader(args):
    beta_contents = Utilities.contentsWithPatternsFromFolder(args.beta_folder, [])
    r = pandas.DataFrame()
    for beta_name in beta_contents:
        logging.info("Processing %s", beta_name)
        beta_path = os.path.join(args.beta_folder, beta_name)
        b = pandas.read_table(beta_path)
        r = pandas.concat([r, b])
    return r

def _gwas_wrapper(gwas):
    logging.info("Processing input gwas")
    return gwas

def build_context(args, gwas):
    logging.info("Loading model from: %s", args.model_db_path)
    model = load_model(args.model_db_path)

    covariance_manager = MatrixManager.load_matrix_manager(args.covariance)

    gwas = _gwas_wrapper(gwas) if gwas is not None else _beta_loader(args)
    context = _build_context(args.model_db_path,args.covariance,model, gwas,covariance_manager)
    return context

def _build_context(model_db_path,covariance,model, gwas,covariance_manager):
    gwas = _prepare_gwas(gwas)
    gwas = _sanitized_gwas(gwas)
    context = OptimizedContext(model_db_path,covariance,gwas, model,covariance_manager)
    return context


def _to_int(d):
    r = d
    try:
        r = int(d)
    except:
        pass
    return r

def format_output(results, context, remove_ens_version):
    results = results.drop("n_snps_in_model",1)
    results = results.drop("n_snps_in_cov",1)
    results = results.drop("n_snps_used",1)
    # Dodge the use of cdf on non finite values
    i = numpy.isfinite(results.zscore)
    results[Constants.PVALUE] = numpy.nan
    results.loc[i, Constants.PVALUE] = 2 * stats.norm.sf(numpy.abs(results.loc[i, Constants.ZSCORE].values))
    return results
