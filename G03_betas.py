#! /usr/bin/env python
"""
This module contains high level code to parse GWAS/GWAMA summary statistics files and parse them into a standard format.
You can use it as stand alone tool to align it to Predictdb Models or jus tconvert the format,
or it can be called from another script to load the data into memory

TODO:
    Implement "streaming read" where the whole file is not read at once, and output as data is read from file.
    (mostly meant to make the unfiltered operation of this script have less memory footprint)

"""
__author__ = 'Donghui'
import gaact
__version__ = gaact.__version__
import logging
import os
import re
import pandas

from timeit import default_timer as timer

from gaact import Constants
from gaact.gwas import GWAS
from gaact.gwas import Utilities as GWASUtilities
from gaact import PredictionModel
from gaact import Utilities
from gaact import Logging
from gaact import Exceptions

def align_data_to_alleles(data, base, left_on, right_on):
    EA, NEA = Constants.EFFECT_ALLELE, Constants.NON_EFFECT_ALLELE
    EA_BASE, NEA_BASE = EA+"_BASE", NEA+"_BASE"
    merged = pandas.merge(data, base, left_on=left_on, right_on=right_on, suffixes=("", "_BASE"))

    alleles_1 = pandas.Series([set(e) for e in zip(merged[EA], merged[NEA])])
    alleles_2 = pandas.Series([set(e) for e in zip(merged[EA_BASE], merged[NEA_BASE])])
    eq = alleles_1 == alleles_2
    merged = merged[eq]

    flipped = merged[EA] != merged[EA_BASE]
    Z = Constants.ZSCORE
    if Z in merged:
        merged.loc[flipped, Z] = - merged.loc[flipped, Z]
    B = Constants.BETA
    if B in merged:
        merged.loc[flipped, B] = - merged.loc[flipped, B]

    merged.loc[flipped, EA] = merged.loc[flipped, EA_BASE]
    merged.loc[flipped, NEA] = merged.loc[flipped, NEA_BASE]

    return merged

def build_betas(args, model, gwas_format, name):
    logging.info("Building beta for %s and %s", name, args.model_db_path if args.model_db_path else "no database")
    load_from = os.path.join(args.gwas_folder, name)
    if model or args.skip_until_header:
        snps = model.snps() if model else None
        snp_column_name = args.snp_column if model else None
        load_from = GWASUtilities.gwas_filtered_source(load_from, snps=snps, snp_column_name=snp_column_name, skip_until_header=args.skip_until_header, separator=args.separator)
    sep = '\s+' if args.separator is None else args.separator
    b = GWAS.load_gwas(load_from, gwas_format, sep=sep, input_pvalue_fix=args.input_pvalue_fix)

    if model is not None:
        PF = PredictionModel.WDBQF
        base = model.weights[[PF.K_RSID, PF.K_EFFECT_ALLELE, PF.K_NON_EFFECT_ALLELE]].drop_duplicates()
        b = align_data_to_alleles(b, base, Constants.SNP, PF.K_RSID)

    b = b.fillna("NA")
    keep = [GWAS.SNP, GWAS.ZSCORE]
    if GWAS.BETA in b: keep.append(GWAS.BETA)
    b = b[keep]
    return b

def validate(args):
    if not args.gwas_folder: raise Exceptions.InvalidArguments("You need to provide an input folder containing GWAS files")

def run(args):
    start = timer()
    validate(args)
    regexp = re.compile(args.gwas_file_pattern) if args.gwas_file_pattern else  None
    names = Utilities.contentsWithRegexpFromFolder(args.gwas_folder, regexp)
    names.sort()

    if len(names) == 0:
        msg = "No GWAS files found on %s with pattern %s" % (args.gwas_folder, args.gwas_file_pattern,)
        raise Exceptions.ReportableException(msg)

    gwas_format = GWASUtilities.gwas_format_from_args(args)
    GWAS.validate_format_basic(gwas_format)
    GWAS.validate_format_for_strict(gwas_format)
    model = PredictionModel.load_model(args.model_db_path) if args.model_db_path else None

    if args.output_folder:
        if not os.path.exists(args.output_folder):
            os.makedirs(args.output_folder)

        for name in names:
            output_path = os.path.join(args.output_folder, name)
            if not ".gz" in output_path:
                output_path += ".gz"
            if os.path.exists(output_path):
                logging.info("%s already exists, delete it if you want it to be done again", output_path)
                continue

            b = build_betas(args, model, gwas_format, name)
            c = "gzip" if ".gz" in name else None
            b.to_csv(output_path, sep="\t", index=False, compression=c)
        end = timer()
        logging.info("Successfully ran GWAS input processing in %s seconds" %(str(end - start)))
    else:
        r = pandas.DataFrame()
        for name in names:
            b = build_betas(args, model, gwas_format, name)
            r = pandas.concat([r,b])
        end = timer()
        logging.info("Successfully parsed input gwas in %s seconds"%(str(end-start)))
        return r

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='M03_betas.py %s: Build betas from GWAS data as expected by MetaXcan.' % (__version__))

    parser.add_argument("--model_db_path",
                        help="Name of model db in data folder. "
                             "If supplied, will filter input GWAS snps that are not present; this script will not produce output if any error is encountered."
                             "If not supplied, will convert the input GWASas found, one line at a atime, until finishing or encountering an error.",
                        default=None)

    parser.add_argument("--gwas_folder",
                        help="name of folder containing GWAS data. All files in the folder are assumed to belong to a single study.",
                        default="data/GWAS")

    parser.add_argument("--gwas_file_pattern",
                        help="Pattern to recognice GWAS files in folders (in case there are extra files and you don't want them selected).",
                        default=None)

    parser.add_argument("--output_folder",
                        help="name of folder to put results in",
                        default=None)

    GWASUtilities.add_gwas_arguments_to_parser(parser)

    parser.add_argument("--separator",
                        help="Character or string separating fields in input file. Defaults to any whitespace.",
                        default=None)

    parser.add_argument("--skip_until_header",
                        help="Some files may be malformed and contain unespecified bytes in the beggining."
                             " Specify this option (string value) to identify a header up to which file contents should be skipped.",
                        default=None)

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    parser.add_argument("--throw",
                        action="store_true",
                        help="Throw exception on error",
                        default=False)

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))
    if args.throw:
        run(args)
    else:
        try:
            run(args)
        except Exceptions.ReportableException as e:
            logging.error("Error:%s", e.msg)
        except Exception as e:
            logging.info("Unexpected error: %s" % str(e))
            exit(1)

