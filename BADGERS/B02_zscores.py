import badgers
__version__ = badgers.__version__

import logging
import os
import numpy
import csv
from timeit import default_timer as timer

from badgers import Logging
from badgers import Utilities
from badgers import Exceptions
from badgers.badgers_help import AssociationCalculation
from badgers.badgers_help import Utilities as badgersUtilities
from badgers.PredictionModel import  WDBQF


def run(args, _gwas=None):
    start = timer()
    logging.info("Started BADGEARS association")
    context = badgersUtilities.build_context(args, _gwas)
    model_snps = context.get_model_snps()
    total_snps = len(model_snps)
    snps_found=set()
    results = []
    logging.info(total_snps)
    trait = args.trait_name
    r, snps = AssociationCalculation.association(trait, context, return_snps=True)
    results.append(r)
    snps_found.update(snps)
    Utilities.ensure_requisite_folders(args.output_file)

    results = AssociationCalculation.dataframe_from_results(zip(*results))
    results = badgersUtilities.format_output(results, context, args.remove_ens_version)
    results.to_csv(args.output_file, index=False, header = False, mode = "a+",quoting=csv.QUOTE_NONNUMERIC)
    end = timer()
    logging.info("Sucessfully processed BADGERS association in %s seconds"%(str(end - start)))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='M04_zscores.py %s: Build ZScores from GWAS data.' % (__version__,))

    parser.add_argument("--model_db_path",
                        help="name of weight db in data folder",
                        default=None)
    parser.add_argument("--input_folder",
                        help="1kg_data_folder",
                        default=None)

    parser.add_argument("--covariance",
                        help="name of file containing covariance data",
                        default=None)

    parser.add_argument("--beta_folder",
                        help="name of folder containing beta data",
                        default=None)

    parser.add_argument("--output_file",
                        help="name of output file",
                        default="results/zscores.csv")

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    parser.add_argument("--remove_ens_version",
                        help="If set, will drop the -version- postfix in gene id.",
                    action="store_true",
                    default=False)

    parser.add_argument("--overwrite",
                        help="If set, will overwrite the results file if it exists.",
                    action="store_true",
                    default=False)

    parser.add_argument("--throw",
                        action="store_true",
                        help="Throw exception on error",
                        default=False)
    parser.add_argument("--trait_name",
                        help="name of trait",
                        default=None)



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
