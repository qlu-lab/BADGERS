#! /usr/bin/env python

import badgers
__version__ = badgers.__version__
import logging
import os

from badgers import Exceptions
from badgers import Logging
from badgers.gwas import Utilities as GWASUtilities

import B01_betas
import B02_zscores_mult


def run(args):
    if not args.model_db_path:
        logging.info("Need to provide a model database file path")
        return
    args.output_folder = None
    g = B01_betas.run(args)
    part_zscore, part_effect,ti = B02_zscores_mult.run(args, g)
    return part_zscore, part_effect,ti

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='BADGERS.py %s:  Will estimate BADGERS results from a set of snp covariance matrices, a model database, and GWAS beta files.' % (__version__))

#weight db model
    parser.add_argument("--model_db_path",
                        help="name of model db in data folder",
                        default=None)

#GWAS betas
    parser.add_argument("--gwas_folder",
                        help="name of folder containing GWAS data. All files in the folder are assumed to belong to a single study.",
                        default="data/GWAS")

    parser.add_argument("--gwas_file_pattern",
                        help="Pattern to recognice GWAS files in folders (in case there are extra files and you don't want them selected).",
                        default=None)

    GWASUtilities.add_gwas_arguments_to_parser(parser)

    parser.add_argument("--separator",
                        help="Character or string separating fields in input file. Defaults to any whitespace.",
                        default=None)

    parser.add_argument("--skip_until_header",
                        help="Some files may be malformed and contain unespecified bytes in the beggining."
                             " Specify this option (string value) to identify a header up to which file contents should be skipped.",
                        default=None)

# ZScore calculation
    parser.add_argument("--covariance",
                        help="name of file containing covariance data",
                        default=None)
                        #default="intermediate/1000GP_Phase3_chr_cov")

    parser.add_argument("--output_file",
                        help="name of output file",
                        default="results/zscores.csv")

    parser.add_argument("--remove_ens_version",
                        help="If set, will drop the -version- postfix in gene id.",
                    action="store_true",
                    default=False)

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    parser.add_argument("--throw",
                        action="store_true",
                        help="Throw exception on error",
                        default=False)

    parser.add_argument("--trait_name",
                        help="name of trait",
                        default=None)
    parser.add_argument("--input_folder",
                        help="1kg_data_folder",
                        default=None)


    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    if args.throw:
        run(args)
    else:
        try:
            run(args)
        except Exceptions.ReportableException as e:
            logging.error(e.msg)
        except Exception as e:
            logging.info("Unexpected error: %s" % str(e))
            exit(1)
