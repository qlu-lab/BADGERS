import B05_database
import re
import os
import shutil

def run(input_file,beta_path,db_path):
    fw = open(input_file,"r")
    with open(input_file, 'r') as fr:
        for text in fr:
            tmp = text.split()
            output = db_path+'/'+tmp[0]+'.db'
            betas = beta_path+'/'+tmp[0]+'.txt'
            B05_database.make_sqlite_db(betas,output,None)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='get last two step done')
    parser.add_argument("--input_file",
                        required=True,
                        help="Path to desired input")
    parser.add_argument("--beta_path",
                        help="path of beta folder",
                        default=None)
    parser.add_argument("--db_path",
                        help="path of db file folder",
                        default=None)

    args = parser.parse_args()
    
    run(args.input_file,args.beta_path,args.db_path)
        
        
            


        
