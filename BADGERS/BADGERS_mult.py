import re
import os
import shutil
import badgers.Formats as Formats
import B03_single_mult as B03_single
import numpy
from scipy import stats


class Args_2(object):
    def __init__(self, args):
        self.model_db_path = args.model_db_path
        self.output_file = args.output_file
        self.gwas_folder = args.gwas_folder
        self.gwas_file_pattern = args.gwas_file_pattern
        self.snp_column = args.snp_column
        self.effect_allele_column = args.effect_allele_column
        self.non_effect_allele_column = args.non_effect_allele_column
        self.zscore_column = args.zscore_column
        self.trait_name = args.trait_name
        self.pvalue_column = args.pvalue_column
        self.beta_column = args.beta_column
        self.separator = None
        self.skip_until_header = None
        self.remove_ens_version = False
        self.verbosity = "10"
        self.throw = False
        self.beta_sign_column = None
        self.chromosome_column = None
        self.position_column = None
        self.freq_column = None
        self.or_column = args.or_column
        self.se_column = args.se_column
        self.input_pvalue_fix = 1e-50
        self.covariance = args.covariance

def run(args):
    if os.path.exists(args.output_file):
        print("%s already exists, move it or delete it if you want it done again", args.output_file)
        return
                
    with open(args.output_file,'wb') as fa:
        fa.write("trait,zscore,effect_size,pvalue\n")
    fa.close
    dest = args.output_file
    head, tail = os.path.split(args.gwas_path)
    args.gwas_folder = head
    args.gwas_file_pattern = tail
    db_path = args.model_db_path
    db_names = []
    trait_names = []
    with open(args.db_name,"r") as fk:
        for lines in fk:
            tmp = lines.split(",")
            db_names.append(tmp[0])
            trait_names.append(lines[lines.index(",")+1:-1])
    covariance_path = args.covariance
    trait_name_list = []
    part_zscore_list = []
    part_effect_list = []
    ti_list = []
    for file in db_names:
        name = file
        trait_number = file.split('.')[0]
        args.covariance = covariance_path+'/'+trait_number+'_cov.txt.gz'
        args.trait_name = str(trait_names[db_names.index(name)])
        trait_name_list.append(args.trait_name)
        args.model_db_path = db_path+'/'+name +'.db'
        args_2 = Args_2(args)
        part_zscore, part_effect,ti = B03_single.run(args_2)
        part_zscore_list.append(part_zscore)
        part_effect_list.append(part_effect)
        ti_list.append(ti)
        print("finish"+name)            
    part_zscore_list = numpy.array(part_zscore_list).transpose()
    part_effect_list = numpy.array(part_effect_list).transpose()
    ti_list = numpy.matrix(ti_list)
    cov_ti = numpy.cov(ti_list)
    cor_ti = numpy.corrcoef(ti_list)
    cov_ti_test = numpy.dot(ti_list,ti_list.transpose())
    cov_ti_inv = numpy.linalg.inv(cov_ti)
    size = len(trait_name_list)
    for trait in trait_name_list:
        I_list = []
        index_now = trait_name_list.index(trait)
        for i in xrange(size):
            if i == index_now:
                I_list.append(float(1))
            else:
                I_list.append(float(0))
        I_list = numpy.array(I_list)
        I_list = numpy.transpose(I_list)
        effect_size = numpy.matmul(numpy.matmul(I_list,cov_ti_inv),part_effect_list)
        zscore = numpy.matmul(numpy.matmul(I_list,cov_ti_inv),part_zscore_list)/numpy.sqrt(cov_ti_inv[index_now][index_now])
        pvalue = 2 * stats.norm.sf(numpy.abs(zscore))
        result = ','.join([str(trait),str(zscore),str(effect_size),str(pvalue)])+"\n"

        with open(dest,'a+') as file:
            file.write(result)
        file.close


    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Build BADGERS.')
#weight db model
    parser.add_argument("--db_name",
                        help="name of traits that going for BADGERS test",
                        default=None)
    
    parser.add_argument("--model_db_path",
                        help="name of model db in data folder",
                        default=None)

#GWAS betas
    parser.add_argument("--gwas_path",
                        help="name of folder containing GWAS data. All files in the folder are assumed to belong to a single study.",
                        default="data/GWAS")


# ZScore calculation

    parser.add_argument("--output_file",
                        help="name of output file",
                        default="results/zscores.csv")


    parser.add_argument("--snp_column",
                        help="SNP_column",
                        default = None)

    parser.add_argument("--effect_allele_column",
                        help="effect_column",
                        default=None)

    parser.add_argument("--non_effect_allele_column",
                        help="non_effect_allele_column",
                        default=None)

    parser.add_argument("--beta_column",
                        help="beta_column",
                        default=None)

    parser.add_argument("--pvalue_column",
                        help="pvalue_column",
                        default=None)

    parser.add_argument("--trait_name",
                        help="name of trait",
                        default="TRAIT")

    parser.add_argument("--gwas_folder",
                        help="name of folder containing GWAS data. All files in the folder are assumed to belong to a single study.",
                        default="data/GWAS")

    parser.add_argument("--gwas_file_pattern",
                        help="Pattern to recognice GWAS files in folders (in case there are extra files and you don't want them selected).",
                        default=None)

    parser.add_argument("--covariance",
                        help="name of file containing covariance data",
                        default=None)
    parser.add_argument("--zscore_column",
                        help="value of zscore",
                        default=None)
    parser.add_argument("--or_column",
                        help="value of or",
                        default=None)
    parser.add_argument("--se_column",
                        help="value of se",
                        default=None)



    args = parser.parse_args()
    
    run(args)

