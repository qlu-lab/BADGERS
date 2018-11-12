
import sqlite3
import os
import Exceptions

class GeneEntry:
    def __init__(self, gene, n_snps):
        self.cho = cho
        self.n_snps = n_snps

class WeightDBEntry:
    def __init__(self, rsid=None, cho=None, weight=None, ref_allele=None, eff_allele=None, pval=None, N=None, cis=None):
        self.rsid = rsid
        self.cho = cho
        self.weight = weight
        self.ref_allele = ref_allele
        self.eff_allele = eff_allele

class WDBQF(object):
    RSID=0
    CHO=1
    WEIGHT=2
    REF_ALLELE=3
    EFF_ALLELE=4


class WeightDB(object):
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
                raise RuntimeError("Weight file doesn't exist")
            self.connection = sqlite3.connect(self.file_name)
            self.cursor = self.connection.cursor()

    def closeDB(self):
        if self.connection:
            self.connection.close()
            self.connection = None
            self.cursor = None

    def weightEntriesFromResults(self, results, result_callback=None):
        weights = []
        for result in results:
            weight = WeightDBEntry(result[WDBQF.RSID],
                                   result[WDBQF.CHO],
                                   result[WDBQF.WEIGHT],
                                   result[WDBQF.REF_ALLELE],
                                   result[WDBQF.EFF_ALLELE])
            weights.append(weight)
            if result_callback:
                result_callback(weight)
        return weights

    def loadFromDB(self, callback=None, cho_key=None):
        self.openDBIfNecessary()

        if cho_key is None:
            results = self.cursor.execute("SELECT rsid, cho, weight, ref_allele, eff_allele FROM weights;")
        else:
            results = self.cursor.execute("SELECT rsid, cho, weight, ref_allele, eff_allele FROM weights where gene = ?;", (gene_key))

        weights = self.weightEntriesFromResults(results, callback)
        return  weights


    def loadGeneNamesFromDB(self):
        self.openDBIfNecessary()
        names = []

        results = self.cursor.execute("SELECT DISTINCT cho FROM weights;")

        for result in results:
            name = result[0]
            names.append(name)

        return names

class WeightDBEntryLogic(object):
    def __init__(self, db_file_name):
        self.weights_by_gene = {}
        self.genes_for_an_rsid = {}
        self.gene_data_for_gene = {}
        self._loadData(db_file_name)

    def anEntryWithRSID(self, rsid):
        entry = None
        if not rsid in self.genes_for_an_rsid:
            return entry

        genes = self.genes_for_an_rsid[rsid]
        gene = genes[0]
        weights = self.weights_by_gene[gene]
        entry = weights[rsid]
        return  entry

    def _loadData(self, db_file_name):
        weights_db = WeightDB(db_file_name)

        class ByNameCallback(object):
            """Helper class to group weights by gene name"""
            def __init__(self, weights_by_gene, genes_for_an_rsid, gene_data_for_gene):
                self.weights_by_gene = weights_by_gene
                self.genes_for_an_rsid = genes_for_an_rsid
                self.gene_data_for_gene = gene_data_for_gene

            def __call__(self, weight):
                if weight.cho in self.weights_by_gene:
                    weights = self.weights_by_gene[weight.cho]
                else:
                    weights = {}
                    self.weights_by_gene[weight.cho] = weights
                weights[weight.rsid]= weight

                if not weight.rsid in self.genes_for_an_rsid:
                    self.genes_for_an_rsid[weight.rsid] = []
                genes = self.genes_for_an_rsid[weight.rsid]

                if not weight.cho in genes:
                    genes.append(weight.cho)


        callback = ByNameCallback(self.weights_by_gene, self.genes_for_an_rsid, self.gene_data_for_gene)
        weights_db.loadFromDB(callback)
