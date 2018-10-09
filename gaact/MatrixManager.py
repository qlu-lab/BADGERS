import pandas
import numpy
import Exceptions
import Formats
class Args(object):
    def __init__(self, weight_db, input_folder):
        self.weight_db = weight_db
        self.covariance_output = None
        self.input_folder = input_folder
        self.verbosity = "10"
        self.input_format = Formats.PrediXcan
        self.min_maf_filter = None
        self.correlation_output = None
        self.max_maf_filter = None

def load_matrix_manager(path):
    d = pandas.read_table(path, sep="\s+")
    m = MatrixManager(d)
    return m


class MatrixManager(object):
    def __init__(self, d):
        _validate(d)
        self.data = _build_data(d)

    def get(self, trait, snps=None, strict=True):
        return _get(self.data, trait, snps, strict)

    def n_snps(self,trait):
        if not trait in self.data:
            return numpy.nan
        snps = self.data[trait]
        snps = _non_na(snps)
        snps = {x[CDTF.RSID1] for x in snps}
        return len(snps)

class CDTF(object):
    TRAIT=0
    RSID1=1
    RSID2=2
    VALUE=3

    K_TRAIT = "TRAIT"
    K_RSID1 = "RSID1"
    K_RSID2 = "RSID2"
    K_VALUE = "VALUE"

def _validate(d):
    processed_traits = set()
    last_trait = None
    traits = d[CDTF.K_TRAIT]
    for g in traits:
        if g != last_trait:
            processed_traits.add(g)
            last_trait = g

   # if numpy.any(d.duplicated()):
   #     msg = "Duplicated SNP entries found"
   #     raise Exceptions.InvalidInputFormat(msg)

def _build_data(d):
    d = d.fillna("NA")
    d = zip(d[CDTF.K_TRAIT].values, d[CDTF.K_RSID1].values, d[CDTF.K_RSID2].values, d[CDTF.K_VALUE].values)
    r = {}
    for t in d:
        trait = 'great'
        if not trait in r:
            r[trait] = []
        r[trait].append(t)
    return r

def _get(d, trait, snps_whitelist=None, strict=True):
    d = d['great']
    if snps_whitelist is not None:
        g, r1, r2, v = zip(*d)
        snps = set(r1)
        snps_whitelist = set(snps_whitelist)
        d = [x for x in d if x[CDTF.RSID1] in snps_whitelist]
    _s_1 = set()
    snps_1 = []
    entries_1 = {}
    for row in d:
        rsid1 = row[CDTF.RSID1]
        rsid2 = row[CDTF.RSID2]
        if rsid1 != rsid2:
            continue
        value = row[CDTF.VALUE]
        if value == "NA":continue
        entries_1[rsid1] = value
        if not rsid1 in _s_1:
            _s_1.add(rsid1)
            snps_1.append(rsid1)
    snps_1 = list(snps_1)
    rows_1 = []
    for snp_i in snps_1:
        rows_1.append(entries_1[snp_i])
    

    return snps_1, rows_1
def _non_na(snps_data):
    return [x for x in snps_data if  x[CDTF.VALUE] != "NA"]
