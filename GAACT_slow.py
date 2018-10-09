import re
import os
import shutil
import gaact.Formats as Formats
import GAACT_single

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
    #if os.path.exists(args.output_file):
       # print("%s already exists, move it or delete it if you want it done again", args.output_file)
        #return
    exist = False
    if os.path.exists(args.output_file):
        exist = True
    if(not exist):
        with open(args.output_file,'wb') as fa:
            fa.write("trait,zscore,effect_size,var_g,pvalue\n")
        fa.close
    wtf = []
    i = 0
    if(exist):
        with open(args.output_file,'r') as fn:
            for line in fn:
                i = i+1
        fn.close

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
    k = 0
    for file in db_names:
        if k<i-1:
            k = k+1
            continue
        name = file
        trait_number = file.split('.')[0]
        args.covariance = covariance_path+'/'+trait_number+'_cov.txt.gz'
        #print(args.covariance)
        args.trait_name = str(trait_names[db_names.index(name)])

        #print(args.trait_name)
        args.model_db_path = db_path+'/'+name
        #print(args.model_db_path)
        args_2 = Args_2(args)
        GAACT_single.run(args_2)
        print("finish"+name)            
    

    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Build GAACT.')
#weight db model
    parser.add_argument("--db_name",
                        help="name of traits that going for GAACT test",
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

