#!/usr/bin/python

import sys, getopt
import os.path
import params
from stitch import stitch

def paras_in_file(paras_file):
    if not os.path.isfile(paras_file):    
        print "ERROR! Parameter file " + paras_file + " doesn't exist!"
        exit()
    argv = []
    for line in open(paras_file):
        line = line.strip()
        cols = line.split()
        argv += cols
    return argv   


def abs_path(path):
    cols = path.split('/')
    if cols[0] == '~':
        user = os.path.expanduser('~') 
        return user + path[1:]
    if cols[0] == '.':
        current = os.path.abspath(os.path.join(os.getcwd(), "."))
        return current + path[1:]
    if cols[0] == '..':
        parent = os.path.abspath(os.path.join(os.getcwd(), ".."))
        return parent + path[2:]
    return path

def main():

    #default input and output
    paras_file = abs_path("./paras.txt")
    fasta_list_file = abs_path("./input_files_list.txt")
    output_dir = abs_path("./output_dir")
    optmap_list_file = abs_path("./optmap_files_list.txt")
   
    #default tools 
    REFALIGNER = "RefAligner" 
    BLASTN = "blastn"
    GLPSOL = "glpsol" 
    DIGEST = abs_path("./fa2cmap_multi.pl")


    #obtaining parameters
    paras_string = "f:i:o:m:t:p:x:y:z:g:a:b:c:d:e:hr"
    opts, args = getopt.getopt(sys.argv[1:], paras_string)
    for op, value in opts:
        if op == "-f":
            paras_file = value
            argv_in_file = paras_in_file(paras_file)
            opts_in_file, args_in_file = getopt.getopt(argv_in_file, paras_string) 
            for fop, fvalue in opts_in_file:
                if fop == "-i":
                    fasta_list_file = abs_path(fvalue)
                elif fop == "-o":
                    output_dir = abs_path(fvalue)
                elif fop == "-m":
                    optmap_list_file = abs_path(fvalue)
                elif fop == "-p":
                    params.num_threads = int(fvalue)
                elif fop == "-x":
                    REFALIGNER = abs_path(fvalue)
                elif fop == "-y":
                    DIGEST = abs_path(fvalue) 
                elif fop == "-g":
                    GLPSOL = abs_path(fvalue)
                elif fop == "-a":
                    params.threshold_1 = float(fvalue)
                elif fop == "-b":
                    params.threshold_2 = float(fvalue)
                elif fop == "-h":
                    params.repeat_naive = False
                elif fop == "-r":
                    params.gap_naive = False

        elif op == "-i":
            fasta_list_file = abs_path(value)
        elif op == "-o":
            output_dir = abs_path(value)
        elif op == "-m":
            optmap_list_file = abs_path(value)
        elif op == "-p":
            params.num_threads = int(value)
        elif op == "-x":
            REFALIGNER = abs_path(value)
        elif op == "-g":
            GLPSOL = abs_path(value)
        elif op == "-y":
            DIGEST = abs_path(value)
        elif op == "-a":
            params.threshold_1 = float(value)
        elif op == "-b":
            params.threshold_2 = float(value)
        elif op == "-h":
            params.repeat_naive = False
        elif op == "-r":
            params.gap_naive = False
    stitch(paras_file, fasta_list_file, output_dir, REFALIGNER,  GLPSOL, DIGEST, optmap_list_file)        

if __name__ == "__main__":
    main()
