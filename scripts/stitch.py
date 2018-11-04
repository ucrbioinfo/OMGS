import os.path
from shutil import copyfile

from merge_inputs import merge_inputs
from chname_fasta import chname_fasta
from fasta_long_seqs import fasta_long_seqs
from core_operations import core_operations
import params

def output_paras(fasta_list_file, output_dir, optmap_list_file, REFALIGNER,  GLPSOL, DIGEST):
    fo = file(output_dir+"/parameters.log", 'w')
    fo.write("fasta list files:\t"+fasta_list_file+"\n")
    fo.write("output directory:\t"+output_dir+"\n")
    fo.write("optmap list files:\t"+optmap_list_file+"\n")
    fo.write("number of threads:\t"+str(params.num_threads)+"\n")
    fo.write("REFALIGNER tool:\t"+REFALIGNER+"\n")
    fo.write("GLPSOL tool:\t"+GLPSOL+"\n")
    fo.write("DIGEST tool:\t"+DIGEST+"\n")
    fo.write("threshold_1:\t"+str(params.threshold_1)+"\n")
    fo.write("threshold_2:\t"+str(params.threshold_2)+"\n")
    fo.write("repeat_naive:\t"+str(params.repeat_naive)+"\n")
    fo.write("gap_naive:\t"+str(params.gap_naive)+"\n")
    fo.close()

def get_optmap_info(optmap_list_file, output_dir):
    fo = file(output_dir+"/optmap_info.log", 'w')
    optmap_list = []
    optmap_type_list = []
    for line in open(optmap_list_file):
        line = line.strip()
        cols = line.split()
        optmap_type = cols[0]
        optmap = cols[1]
        fo.write(optmap_type+"\t"+optmap+"\n")
        optmap_list.append(optmap)
        optmap_type_list.append(optmap_type)
    return optmap_list, optmap_type_list

def stitch(paras_file, fasta_list_file, output_dir, REFALIGNER, GLPSOL, DIGEST, optmap_list_file):


    #intermediate files unrelated to each iteration
    input_merged_file = output_dir + "/input_merged.fasta" 
    input_merged_file_chname = output_dir + "/input_merged_chname.fasta"
    input_merged_file_chname_long = output_dir + "/input_merged_chname_long.fasta"
    sub_dir = output_dir + "/core"
    

    #check existence of output_dir
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    #pre-operations
    output_paras(fasta_list_file, output_dir, optmap_list_file, REFALIGNER,  GLPSOL, DIGEST) 
    merge_inputs(output_dir, fasta_list_file, input_merged_file)
    chname_fasta(input_merged_file, input_merged_file_chname)
    fasta_long_seqs(input_merged_file_chname, input_merged_file_chname_long) 
    optmap_list, optmap_type_list = get_optmap_info(optmap_list_file, output_dir)

    refaligner_dir = sub_dir + "/refaligner"
    contigs_dir = sub_dir + "/grouped_contigs"
    listfiles_dir = sub_dir + "/grouped_listfiles"

    input_fasta_file = sub_dir + "/input_contigs.fasta"
    mtp_fasta_file = sub_dir + "/input_contigs_mtp.fasta"
    mtp_fasta_file_chname = sub_dir + "/input_contigs_mtp_chname.fasta"
    unaligned_fasta_file = sub_dir + "/input_contigs_unaligned.fasta"
    lowconf_fasta_file = sub_dir + "/input_contigs_lowconf.fasta"

    refaligner_prefix = "input_contigs"

        
    #pre-processing
    if not os.path.isdir(sub_dir):
        os.makedirs(sub_dir)
    copyfile(input_merged_file_chname_long, input_fasta_file) 
        
    #core operations
    core_operations(REFALIGNER, DIGEST, GLPSOL, sub_dir, refaligner_dir, input_fasta_file, optmap_list, optmap_type_list, refaligner_prefix, mtp_fasta_file, mtp_fasta_file_chname, unaligned_fasta_file, lowconf_fasta_file)
         

