#!/usr/bin/env python

import sqlite3
import sys
from timeit import default_timer as timer

class DB:
    def __init__(self, db_file):
        '''Opens connection to database'''
        self.connection = sqlite3.connect(db_file)
        self.c = self.connection.cursor()
        self.weights_fields = {"rsid", "cho", "weight", "ref_allele", "eff_allele"}
        self.extra_fields = {"gene", "genename", "pred.perf.R2", "n.snps.in.model", "pred.perf.pval", "pred.perf.qval"}
        self.construction_fields = {"chr", "cv.seed"}
        self.meta_fields = {"n.samples"}

    def add_weights(self, betas):
        # Drop weights table if it already exists.
        self.c.execute("DROP INDEX IF EXISTS weights_rsid_gene")
        self.c.execute("DROP INDEX IF EXISTS weights_gene")
        self.c.execute("DROP INDEX IF EXISTS weights_rsid")
        self.c.execute("DROP TABLE IF EXISTS weights")
        # Create new weights table.
        self.c.execute('CREATE TABLE weights ' +
            '(rsid TEXT, cho TEXT, weight DOUBLE, ref_allele CHARACTER, ' + 
            'eff_allele CHARACTER)')
        self.c.execute("CREATE INDEX weights_rsid ON weights (rsid)")
        self.c.execute("CREATE INDEX weights_gene ON weights (cho)")
        self.c.execute("CREATE INDEX weights_rsid_gene ON weights (rsid, cho)")
        # Insert data.
        for row in data_rows(betas, self.weights_fields):
            self.insert_weights_row(row)

    def insert_weights_row(self, row):
        # Insert a row into the weights table. Variable row is a dict.
        self.c.execute("INSERT INTO weights VALUES (?, ?, ?, ?, ?)",
            (row['rsid'], row['cho'], row['weight'], row['ref_allele'],row['eff_allele']))

    def add_extra(self, results):
        # Drop extra table if it already exists.
        self.c.execute("DROP INDEX IF EXISTS extra_gene")
        self.c.execute("DROP TABLE IF EXISTS extra")
        # Create new extra table with index on gene.
        self.c.execute("CREATE TABLE extra (gene TEXT, genename TEXT, " + 
            "`pred.perf.R2` DOUBLE, `n.snps.in.model` INTEGER, `pred.perf.pval` DOUBLE, `pred.perf.qval` DOUBLE)")
        self.c.execute("CREATE INDEX extra_gene ON extra (gene)")
        # Insert data.
        for row in data_rows(results, self.extra_fields):
            self.insert_extra_row(row)

    def insert_extra_row(self, row):
        self.c.execute("INSERT INTO extra VALUES (?, ?, ?, ?, ?, ?)",
            (row['gene'],row['genename'],row['pred.perf.R2'],row['n.snps.in.model'],row['pred.perf.pval'], row['pred.perf.qval']))

    def add_construction(self, construction):
        # Drop construction table if it already exists.
        self.c.execute("DROP TABLE IF EXISTS construction")
        # Create new construction table.
        self.c.execute("CREATE TABLE construction (chr INTEGER, " + 
            "`cv.seed` INTEGER)")
        # Insert data.
        for row in data_rows(construction, self.construction_fields):
            self.insert_construction_row(row)


    def insert_construction_row(self, row):
        self.c.execute("INSERT INTO construction VALUES (?, ?)",
            (row['chr'], row['cv.seed']))

    def add_meta(self, meta):
        # Drop meta_data table if it already exists.
        self.c.execute("DROP TABLE IF EXISTS sample_info")
        # Create new meta_data table.
        self.c.execute("CREATE TABLE sample_info (`n.samples` INTEGER)")
        # Insert data.
        for row in data_rows(meta, self.meta_fields):
            self.insert_meta_row(row)

    def insert_meta_row(self, row):
        self.c.execute("INSERT INTO sample_info VALUES(?)", (row['n.samples'],))

    def close(self):
        self.connection.commit()
        self.connection.close()

class Intron_DB(DB):

    def __init__(self, db_file):
        DB.__init__(self, db_file)
        self.construction_fields.add("pos_mod_100")

    def add_construction(self, construction):
        # Drop construction table if it already exists.
        self.c.execute("DROP TABLE IF EXISTS construction")
        # Create new construction table.
        self.c.execute("CREATE TABLE construction (chr INTEGER, " + 
            "`n.genes` INTEGER, `cv.seed` INTEGER, pos_mod_100 INTEGER)")
        # Insert data.
        for row in data_rows(construction, self.construction_fields):
            self.insert_construction_row(row)


    def insert_construction_row(self, row):
        self.c.execute("INSERT INTO construction VALUES (?, ?, ?, ?)",
            (row['chr'], row['n_genes'], row['seed_for_cv'], row['pos_mod_100']))


def upconvert(x):
    '''Helper function for data_rows().'''
    for f in (int, float):
        try:
            return f(x)
        except ValueError:
            pass
    return x


def data_rows(source_file, expected_fields):
    '''
    Iterate over data rows in the source file, labeling fields and converting
    formats as required.

    source_file - the name of the tab-delimited file to be read from.
    expected_fields - A set of the fields expected to be found the source file.
    '''
    header = None
    with open(source_file, 'r') as src:
        for k, line in enumerate(src):
            if k == 0:
                if not expected_fields == set(line.strip().split()):
                    raise RuntimeError("Invalid header in " + source_file)
                header = line.strip().split()
            else:
                yield dict(zip(header, map(upconvert, line.strip().split())))


def make_sqlite_db(betas, output, intron):
    start = timer() 
    print("Generating %s..." % output)
    db = Intron_DB(output) if intron else DB(output)
    print("\tCreating weights table...")
    db.add_weights(betas)
    print("\tCreating extra table...")
    db.close()
    print("\tDone.")
    end = timer()
    print("Sucessfully processed weight db in %s seconds"%(str(end - start)))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Create a PredictDB database from input files.')
    parser.add_argument(
        "--output",
        required=True,
        help="Path to desired output")
    parser.add_argument(
        "--betas",
        required=True,
        help="Path to input file with beta values.")
    parser.add_argument(
        "--intron",
        action="store_true",
        help="Include if this was intron model")

    args = parser.parse_args()

    make_sqlite_db(args.betas, args.output, args.intron)
