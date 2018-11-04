

# OMGS Software



# DESCRIPTION


OMGS is a fast genome scaffolding tool which takes advantage of one or multiple Bionano optical maps to accurately generate scaffolds. Instead of alternatively using single optical maps, OMGS uses multiple optical maps at the same time and takes advantage of the redundance contained in multiple maps to generate the ”optimal” scaffolds which make the smartest tradeoff between contiguity and correctness. Extensive experimental results demonstrate that our tool can outperform existing method in terms of both correctness and contiguity on multiple Bionano optical maps, and it can also generate comparable or better scaffolds using single optical map.
Until now, Novo&Stitch can only run on Unix/Linux systems.  



# DEPENDENCY


1. python   
The majoy part of OMGS is written in Python, so Python has to be installed. 
Python2.7(or above) is suggested.  

2. perl   
In OMGS, we use one perl script "fa2cmap_multi.pl" of a scaffolding tool Irys-scaffolding.
So perl has to be installed. 

3. GLPK  
GLPK is a tool for solving Linear Programming (LP) model. Since the false alignments removal module and Minimum Disagreement Correlation Clustering problem solvor in OMGS both uses a standalone tool of GLPK called "glpsol" to solve LP model, GLPK has to be installed. 
GLPK can be found from https://www.gnu.org/software/glpk/#TOCdownloading.

4. RefAligner   
RefAligner is a tool developed by BioNano company to align contigs to optical maps. It is called by OMGS, so it has to be installed. 
RefAligner can be found from http://www.bnxinstall.com/RefalignerAssembler/Linux/ 

5. fa2cmap_multi.pl   
In OMGS, we use one perl script "fa2cmap_multi.pl" of a scaffolding tool Irys-scaffolding. 
The users need to download this script from https://github.com/i5K-KINBRE-script-share/Irys-scaffolding/blob/e8e8f177dce2bf59421bd00c517ab7dc683e25d4/KSU_bioinfo_lab/assemble_XeonPhi/third-party/fa2cmap_multi.pl
   

6. python packages

numpy, scipy, joblib



# PARAMETERS


There are two kinds of parameters, functional parameters and performance related parameters. 
Functional parameters are the ones that users may have to set by themselves. For performance related parameters, most users may not be able to understand the meaning of them, unless they are very familiar with every step of Novo&Stitch algorithm. Therefore, we offer two sets of suggested parameters for users in "performance related parameters" section. 

1. functional parameters  

-f: specify the parameter file which lists all the input paramerters.
Novo$Stitch offers two ways to input parameters for users as follow:

(1) first way
e.g. 

$cd ./OMGS

$python ./scripts/main.py -x /home/stelo/BIONANO_in_progress/tools/RefAligner -y /home/weihua/Novo_Chimeric/tools/fa2cmap_multi.pl -i /home/weihua/cowpea/fastas_cowpea_real_2_quiver91x_high_outcov100.txt -o /home/weihua/cowpea/novo_scaffording_BspQIBssSI_q0_new -m /home/weihua/cowpea/optmaps_for_cowpea.txt -p 32 -a 10 -r 

(2) second way  
e.g.

$cd ./OMGS

$python ./scripts/main.py -f ~/cowpea/parameters.txt

Then in parameters.txt, the parameters are listed line by line as follows:

-x /home/stelo/BIONANO_in_progress/tools/RefAligner 

-y /home/weihua/Novo_Chimeric/tools/fa2cmap_multi.pl 

-i /home/weihua/cowpea/fastas_cowpea_real_2_quiver91x_high_outcov100.txt 

-o /home/weihua/cowpea/novo_scaffording_BspQIBssSI_q0_new 

-m /home/weihua/cowpea/optmaps_for_cowpea.txt 

-p 32 

-a 10 

-r 

For all the parameters about address, absolute paths are preferred. If you use relative paths, '.' and '..' always represent './OMGS' directory and its parent directory respectively, and '~' represents current user's repository.  


-i: specify the fasta list file which lists the address of input fasta files line by line.  
e.g. fasta_list.txt 

/home/stelo/phytophthora/seven_2/seven_2.fasta  

/home/weihua/canu.fasta  

./falcon.fasta  

This parameter is required to be specified by users.   

For all the parameters about address, absolute paths are preferred. If you use relative paths, '.' and '..' always represent './OMGS' directory and its parent directory respectively, and '~' represents current user's repository.  

-o: specify the output dirctory which will contains all of the output files and intermediate files generated. This parameter is required to be specified by users.   

-m: specify the optical map file. This parameter is required to be specified by users.   

e.g. optmaps_for_cowpea.txt

BspQI   /home/stelo/BIONANO_in_progress/vu_162_180K.cmap

BssSI   /home/stelo/BIONANO_in_progress/vu_bsss1_102.cmap

This parameter is required to be specified by users.   

For all the parameters about address, absolute paths are preferred. If you use relative paths, '.' and '..' always represent './OMGS' directory and its parent directory respectively, and '~' represents current user's repository.  


-p: specify the number of threads. The default value is "32".  

-x: specify the address of Refaligner executable program on your machine. e.g. /home/weihua/tools/RefAligner
If you can run RefAligner in any directory of your machine by command like "RefAligner [parameters]", then you can put "RefAligner" as this parameter after "-x". Anyway, this parameter should be the first term of the command your input when you run RefAligner on your machine. The default value is "RefAligner".   

-y: specify the address of fa2cmap_multi.pl file. e.g. /home/weihua/OMGS/tools/fa2cmap_multi.pl

-g: specify the address of glpsol executable prgram on your machine. e.g. /home/weihua/tools/glpsol  
You should be able to find it after installing GLPK.  
Same as Refaligner, you may also be able to use just "glpsol". The default value is "glpsol".  


2. performance related parameters

-a: specify threshold_1 which is used in pre-processing. Only the alignments whose confidences from Refaligner are larger than threshold_1 are used. The default value is "15".  

-b: specify threshold_2 which is used in false alignments removal. When building hypergraph, if the diffence between the two distance between two contigs' middle point coordinates on two optical map molecules is larger than threshold_2 of larger contig's length, one hyperedge is built on the four alignments. The default value is "0.2".  

-h: With "-h", OMGS uses statistical test to recognize the repetitive regions in optical moleculs rather than naive method. Statistical test might be able to get more accurate result but is much slower. 

-r: With "-r", OMGS uses statistical model to estimate gaps rather than naive method. Statistical model might be able to get more accurate result but is a little bit slower.



# OUTPUT FILES

There are two kinds of output files, stitched contigs file and log files. Stitched contigs file contains the final stitched contigs, while the log files contain the information which shows the whole process of running. All of the output files are stored in output directory specified by users. There are some other intermediate files inside, which could be ignored by users.
 
1. stitched contigs file  

final_stitched.fasta: It's a fasta file which contains the stitched contigs.  

2. log files and some intermediate files  

(1) in main output directory  

parameters.log: shows the values of all parameters used. 

input.log: shows all the input fasta files  

(2) in each iteration_* directory

alms_0_initial.log: all alignments 

alms_1_removed_lowconf.log: alignments after removing low quilty alignments

alms_2_removed_false_alms.log: alignments after removing false alignments

alms_3_removed_contained_locally.log: alignments after removing contained alignments locally

alms_4_unified.log: alignments after coordinates unification

alms_5_removed_contained_globally.log: alignments after removing contained alignments globally

alms_6_mtp.log: alignments of mtp contigs

ugraph_1.log: first association graph

ugraph_2.log: second association graph

forest_1.log: msf of first association graph

forest_2.log: msf of second association graph

shifts.log: shifts information for coordinates unification

alm_graph_edges.txt: edges of alignment graph

alm_graph_vertices.txt: vertices of alignment graph

cover.log: optimal vertex cover got for false alignments removal

mtp_nodes.log: mtp contigs

stitch_info.log: stitch information

confidences.txt: the confidences of the contigs mapped to stitched contigs

cancel_id_list.txt: the information of cancelled stitch

input_contigs.fasta: input contigs for this iteration

input_contigs_mtp.fasta: mtp contigs for this iteration

adjusted_contigs_total.fasta: ouput stitched contigs for this iteration


