import os
import sqlite3
import pandas

import Exceptions

class WDBQF(object):
    "Weight DB weight Query Format"
    RSID=0
    TRAIT=1
    WEIGHT=2
    REF_ALLELE=3
    EFF_ALLELE=4

    K_RSID="rsid"
    K_TRAIT="cho"
    K_WEIGHT="weight"
    K_EFFECT_ALLELE="effect_allele"
    K_NON_EFFECT_ALLELE="non_effect_allele"

    ORDER = [(K_RSID, RSID), (K_TRAIT,TRAIT), (K_WEIGHT,WEIGHT), (K_EFFECT_ALLELE, EFF_ALLELE),(K_NON_EFFECT_ALLELE, REF_ALLELE)]


def snps_in_db(path):
    t = ModelDB(path)
    w = t.load_weights()
    s = set(w[WDBQF.RSID])
    return s

class ModelDB(object):
    def __init__(self, file_name , create_if_absent=False):
        self.connection = None
        self.cursor = None
        self.file_name = file_name
        self.create_if_absent = create_if_absent

    def __del__(self):
        self.closeDB()

    def openDBIfNecessary(self):
        if not self.connection:
            if not self.create_if_absent and not os.path.exists(self.file_name):
                raise Exceptions.BadFilename(self.file_name)
            self.connection = sqlite3.connect(self.file_name)
            self.cursor = self.connection.cursor()

    def closeDB(self):
        if self.connection:
            self.connection.close()
            self.connection = None
            self.cursor = None

    def load_weights(self, trait_key=None):
        self.openDBIfNecessary()

        params = []
        query, params = query_helper("SELECT rsid, cho, weight, ref_allele, eff_allele FROM weights", trait_key)
        while True:
            try:
                results = self.cursor.execute(query, params)
            except sqlite3.OperationalError as e:
                print("fc")
                continue
                raise Exceptions.ReportableException("Could not read input tissue database. Please try updating the tissue model files.")
            except Exception as e:
                raise e
            break
        weights = zip(*results)
        return  weights

    
def query_helper(query, trait_key=None):
    params = []
    if trait_key:
        query += " WHERE trait = ?"
        params.append(trait_key)
    query += ";"
    return query, tuple(params)

class Model(object):
    def __init__(self, weights):
        self.weights = weights

    def snps(self):
        snps = self.weights.rsid.values
        return set(snps)

def dataframe_from_weight_data(w):
    weights = pandas.DataFrame({key: w[order] for key, order in WDBQF.ORDER})
    weights = weights[[key for key,order in WDBQF.ORDER]]
    return weights


def load_model(path):
    db = ModelDB(path)
    weights = db.load_weights()
    weights = dataframe_from_weight_data(weights)
    model = Model(weights)
    return model
