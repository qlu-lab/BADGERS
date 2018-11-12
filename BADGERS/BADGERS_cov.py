import B04_covariance
import re
import os
import shutil
import badgers.Formats as Formats
import badgers

class Args(object):
    def __init__(self, weight_db, input_folder, covariance_output):
        self.weight_db = weight_db
        self.input_folder = input_folder
        self.covariance_output = covariance_output
        self.verbosity = "10"
        self.input_format = Formats.PrediXcan
        self.min_maf_filter = None
        self.correlation_output = None
        self.max_maf_filter = None
        self.max_snps_in_gene = None
                                    

        
def run(input_file,model_db_path,input_folder,covariance):
    fw = open(input_file,"r")
    for text in fw:
        print(text)
        text = text.split()[0]
        args = Args(model_db_path+"/"+text+".db",input_folder,
                    covariance +"/"+text+"_cov.txt.gz")
        work = B04_covariance.ProcessWeightDB(args)
        work.run()
        print("finish_cov")
    fw.close()
    

    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='get last two step done')
    parser.add_argument(
        "--input_file",
        required=True,
        help="Path to desired input")
    parser.add_argument("--model_db_path",
                        help="name of model db in data folder",
                        default=None)
    parser.add_argument("--input_folder",
                        help="location of 1kg genome data",
                        default=None)
    parser.add_argument("--covariance",
                        help="name of file containing covariance data",
                        default=None)

    args = parser.parse_args()
    
    run(args.input_file,args.model_db_path,args.input_folder,args.covariance)

