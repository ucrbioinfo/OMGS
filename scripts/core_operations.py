import os.path
from shutil import copyfile

from core import core
from scaffolding import scaffolding
from flip_xmap import flip_xmap

def core_operations(REFALIGNER, DIGEST, GLPSOL, sub_dir, refaligner_dir, input_fasta_file, optmap_list, optmap_type_list, refaligner_prefix, mtp_fasta_file, mtp_fasta_file_chname, unaligned_fasta_file, lowconf_fasta_file):
    fasta_file_in_refaligner = refaligner_dir + "/" + refaligner_prefix + ".fasta"

    if not os.path.isdir(refaligner_dir):
        os.makedirs(refaligner_dir)
    copyfile(input_fasta_file, fasta_file_in_refaligner)

    out_list = []
    out2_list = []
    for i in range(len(optmap_list)):
        optmap = optmap_list[i]
        optmap_type = optmap_type_list[i]
        silicomap = refaligner_dir + "/" + refaligner_prefix + "_" + optmap_type + ".cmap"
        out = refaligner_dir + "/" + refaligner_prefix + "_" + optmap_type + "_BNG_VS_seq"
       
        out2 = refaligner_dir + "/" + refaligner_prefix + "_" + optmap_type

        command = "perl " + DIGEST + " -v -i " + fasta_file_in_refaligner + " -e " + optmap_type + " -m 0 -M 0"
        os.system(command)
       
        command = REFALIGNER + " -i " + optmap + " -ref " + silicomap + " -o " + out + " -res 2.9 -FP 1.2 -FN 0.15 -sf 0.10 -sd 0.15 -extend 1 -outlier 1e-4 -endoutlier 1e-2 -deltaX 12 -deltaY 12 -xmapchim 14 -T 1e-8 -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 50 -mres 1e-3 -insertThreads 4 -nosplit 2 -biaswt 0 -indel -rres 1.2 -f -maxmem 256 "
        os.system(command)
        flip_xmap(out)

        out_list.append(out)
        out2_list.append(out2)    
   
    core(optmap_list, input_fasta_file, out_list, out2_list, sub_dir, GLPSOL) 
    orders_file = sub_dir + "/orders.txt"
    gaps_file = sub_dir + "/gaps.txt"
    ori_file = sub_dir + "/contig_oris.txt"
    scaffolding(out2_list, input_fasta_file, orders_file, gaps_file, ori_file, sub_dir)

