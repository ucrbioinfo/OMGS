#!usr/bin/python
import csv
import sys

from sets import Set

from collections import defaultdict

import collections
import math

from alignment import Alignment 
from get_mst import get_mst
from merge_order_graphs import merge_order_graphs
from false_alms import false_alms
from repeats import get_repeats, get_var
from full_order import check_full_order, get_candiedges
import full_order
from gaps import gaps
from minFAS import minFAS
import params

def copy_alms(old_alms, removed):
    current_alms = {}
    for ref in old_alms:
        current_alms[ref] = {}
        for qry in old_alms[ref]:
            x = old_alms[ref][qry]
            if removed[ref, qry] == False:
                current_alms[ref][qry] = x
    return current_alms

def print_alms(alms):
    for ref in alms:
        for qry in alms[ref]:
            x = alms[ref][qry]
            print "alignment:", (ref, qry), x.start, x.end

def count_alms(alms):
    num_alms = 0
    for ref in alms:
        num_alms += len(alms[ref])
    return num_alms

def output_alms(alms, filename):
    fo = file(filename, 'w')
    fo.write("ref_id\tqry_id\torientation\tconfidence\tstart\tend\n")
    for ref in alms:
        for qry in alms[ref]:
            x = alms[ref][qry]
            fo.write(str(ref)+"\t"+str(qry)+"\t"+str(x.orientation)+"\t"+str(x.confidence)+"\t"+str(x.start)+"\t"+str(x.end)+"\n")
    fo.close()

def different_contigs(alms_1, alms_2):
    contigs_1 = set([])
    for ref in alms_1:
        for qry in alms_1[ref]:
            contigs_1.add(qry)
    contigs_2 = set([])
    for ref in alms_2:
        for qry in alms_2[ref]:
            contigs_2.add(qry)
    diff = set([])
    for c in contigs_1:
        if c not in contigs_2:
            diff.add(c)
    return diff

def output_contigs(contigs, filename):
    fo = file(filename, 'w')
    for c in contigs:
        fo.write(str(c)+"\n")
    fo.close()

def output_forest(forest, vertex_orientations, filename):
    fo = file(filename, 'w')
    for root in forest:
        fo.write("root:"+str(root)+"\n")
        oppo_tree = {}
        for v in forest[root]:
            for u in forest[root][v]:
                oppo_tree[u] = v
        oppo_tree[root] = -1

        for v in forest[root]:
            fo.write("node:"+str(v)+"\t")
            fo.write(str(vertex_orientations[root][v])+"\t")
            fo.write("parent:"+str(oppo_tree[v])+"\t")
            fo.write("children: ")
            for u in forest[root][v]:
                fo.write(str(u)+" ")
            fo.write("\n")
        fo.write("\n")
    fo.close()

def search_in_sorted_list(L, v):
    smaller_num = 0
    larger_num = 0
    for i in range(len(L)):
        if L[i] < v:
            smaller_num += 1
        elif L[i] > v:
            larger_num += 1
    return (smaller_num, larger_num)

def markers_in_qry_left_overhang(qry_markers, x):
    qry = x.qry
    pos = x.qrystartpos
    ori = x.orientation
    pos_list = qry_markers[qry]
    (smaller_num, larger_num) = search_in_sorted_list(pos_list, pos)
    if ori == '+':
        return smaller_num
    else:
        return larger_num

def markers_in_qry_right_overhang(qry_markers, x):
    qry = x.qry
    pos = x.qryendpos
    ori = x.orientation
    pos_list = qry_markers[qry]
    (smaller_num, larger_num) = search_in_sorted_list(pos_list, pos)

    if ori == '+':
        return larger_num
    else:
        return smaller_num

def output_DAGs(DAGs, filename):
    fo = file(filename, 'w')
    for i in range(0, len(DAGs)):
        DAG = DAGs[i]
        oppo_DAG = {}
        for v in DAG:
            oppo_DAG[v] = set([])
        for v in DAG:
            for u in DAG[v]:
                oppo_DAG[u].add(v)
        fo.write("source: ")
        for v in oppo_DAG:
            if oppo_DAG[v] == set([]):
                fo.write(str(v)+" ")
        fo.write("\n")
        fo.write("sink: ")
        for v in DAG:
            if DAG[v] == set([]):
                fo.write(str(v)+" ")
        fo.write("\n")
        
        for v in DAG:
            fo.write("node:"+str(v)+"\t")
            fo.write("incoming: ")
            for u in oppo_DAG[v]:
                fo.write(str(u)+" ")
            fo.write("\t")
            fo.write("outgoing: ")
            for u in DAG[v]:
                fo.write(str(u)+" ")
            fo.write("\n")
        fo.write("\n")
    fo.close()


def topolgical_sort(graph_unsorted):
    """
    Repeatedly go through all of the nodes in the graph, moving each of
    the nodes that has all its edges resolved, onto a sequence that
    forms our sorted graph. A node has all of its edges resolved and
    can be moved once all the nodes its edges point to, have been moved
    from the unsorted graph onto the sorted one.
    """

    # This is the list we'll return, that stores each node/edges pair
    # in topological order.
    graph_sorted = []

    # Convert the unsorted graph into a hash table. This gives us
    # constant-time lookup for checking if edges are unresolved, and
    # for removing nodes from the unsorted graph.
    graph_unsorted = dict(graph_unsorted)

    # Run until the unsorted graph is empty.
    while graph_unsorted:

        # Go through each of the node/edges pairs in the unsorted
        # graph. If a set of edges doesn't contain any nodes that
        # haven't been resolved, that is, that are still in the
        # unsorted graph, remove the pair from the unsorted graph,
        # and append it to the sorted graph. Note here that by using
        # using the items() method for iterating, a copy of the
        # unsorted graph is used, allowing us to modify the unsorted
        # graph as we move through it. We also keep a flag for
        # checking that that graph is acyclic, which is true if any
        # nodes are resolved during each pass through the graph. If
        # not, we need to bail out as the graph therefore can't be
        # sorted.
        acyclic = False
        for node, edges in list(graph_unsorted.items()):
            for edge in edges:
                if edge in graph_unsorted:
                    break
            else:
                acyclic = True
                del graph_unsorted[node]
                graph_sorted.append((node, edges))

        if not acyclic:
            # Uh oh, we've passed through all the unsorted nodes and
            # weren't able to resolve any of them, which means there
            # are nodes with cyclic edges that will never be resolved,
            # so we bail out with an error.
            raise RuntimeError("A cyclic dependency occurred")

    return graph_sorted

def pre_process(optmap_i, optmap_file, myfile, myfile2, output_dir, min_confidence):
    header_lines = 10
    header = []
    minrefoverhang = 50000
    minqryoverhang = 50000

    all_alms = {} # stores all the Alignments for all groups, all_groups[ref] should contain molecule ref
    qualify_alms = {} # only keep one alignment(the one with highest confidence) for each contig in one molecule
    removed = {} # removed[ref,qry] == True means alignment for (ref, qry) is already removed

    # collecting alignments and store in all_groups
    print '---------------read .xmap file-------------------'
    with open(myfile+'_flip.xmap', 'rb') as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        for i in range(header_lines): # 10 lines of header
            header.append(csvreader.next()) # save them
        # read the first non-header line
        while True:
            try:
                row = csvreader.next()
                x = Alignment(int(row[1]),int(row[2]),float(row[3]),float(row[4]),float(row[5]),
                                  float(row[6]),row[7],float(row[8]), row[9],float(row[10]),
                                  float(row[11]),int(row[12]),row[13])
                if x.ref not in all_alms:
                    all_alms[x.ref] = [x]
                else:
                    all_alms[x.ref].append(x)
            except StopIteration:
                break
    num_all_alms = 0
    for ref in all_alms:
        num_all_alms += len(all_alms[ref])
    print "In total, the number of alignments collected is ", num_all_alms
    
    # only keep one alignment(the one with highest confidence) for each contig in one molecule    
    for ref in all_alms:
        group = all_alms[ref]
        qry_bestx = {}
        for x in group:
            if x.qry not in qry_bestx:
                qry_bestx[x.qry] = x
            else:
                if x.confidence > qry_bestx[x.qry].confidence:
                    qry_bestx[x.qry] = x

        qualify_alms[ref] = {}
        for qry in qry_bestx:
            qualify_alms[ref][qry] = qry_bestx[qry]

    num_qualify_alms = 0
    for ref in qualify_alms:
        num_qualify_alms += len(qualify_alms[ref])
    # initialize removed array
    for ref in qualify_alms:
        for qry in qualify_alms[ref]:
            removed[ref,qry] = False
    current_alms = copy_alms(qualify_alms, removed)
    output_alms(current_alms, output_dir+"/opt_"+str(optmap_i)+"_alms_0_initial.log")
    print "In total, the number of alignments in qualify_alms is ", num_qualify_alms



    # remove low confidence alignments
    print '---------------Remove low quality alignments---------------'
    for ref in qualify_alms:
        for qry in qualify_alms[ref]:
            x = qualify_alms[ref][qry]
            if x.confidence < min_confidence:
                removed[ref, qry] = True
                print 'alignment (', ref, ',', qry, ') is low quality and removed'
    num_alms = 0
    for ref in qualify_alms:
        for qry in qualify_alms[ref]:
            if removed[ref, qry] == False:
                num_alms += 1
    current_alms = copy_alms(qualify_alms, removed)
    output_alms(current_alms, output_dir+"/opt_"+str(optmap_i)+"_alms_1_removed_low_conf.log")
    print "After removing low confidence alignments, the number of alignments is ", num_alms
    print '---------------End---------------'

    # read optical map
    optmap = {}
    with open(optmap_file) as f_map:
        for line in f_map:
            line = line.strip()
            if line[0] == '#':
                continue 
            cols = line.split('\t')
            CMapId = int(cols[0])
            LabelChannel = cols[4] 
            Position = float(cols[5])
            
            if CMapId not in optmap:
                optmap[CMapId] = []
            if LabelChannel == "1":
                optmap[CMapId].append(Position)
    for CMapId in optmap:
        optmap[CMapId].sort()

   
    print '---------------scaling-------------------'
    # calculating scaling
    qry_len = {}
    with open(myfile2+'_key.txt') as f_key:
        for i in range(0, 4): # 4 header lines
            f_key.readline()  
        for line in f_key:
            line = line.strip()
            cols = line.split('\t')
            qry_id = int(cols[0])
            seq_len = int(cols[2])
            qry_len[qry_id] = seq_len
    scaling = 0
    num = 0
    with open(myfile+'_r.cmap') as f_q:
        for i in range(0, 11): # 11 header lines
            f_q.readline()
        for line in f_q:
            line = line.strip()
            cols = line.split('\t')
            qry_id = int(cols[0])
            appr_len = float(cols[1])
            seq_len = qry_len[qry_id]
            scaling += appr_len/seq_len
            num += 1
    scaling /= num # scaling=1.02258059775
    scaling = 1.0
    # use scaling to adjsut coordinates of alignments
    for ref in qualify_alms:
        for qry in qualify_alms[ref]:
            x = qualify_alms[ref][qry]
            x.qrystartpos /= scaling
            x.qryendpos /= scaling
            x.qrylen /= scaling
            x.refstartpos /= scaling
            x.refendpos /= scaling
            x.reflen /= scaling          

    # use scaling to adjsut coordinates of optial map
    for ref in optmap:
        for i in range(0, len(optmap[ref])):
            optmap[ref][i] /= scaling

    print '---------------END-------------------'

    # find the reference-based coordinates for each contig
    for ref in qualify_alms:
        for qry in qualify_alms[ref]:
            x = qualify_alms[ref][qry]
            if (x.orientation == '+'):
                x.qry_left_overlen = x.qrystartpos
                x.qry_right_overlen = x.qrylen - x.qryendpos
            else:
                x.qry_left_overlen = x.qrylen - x.qrystartpos
                x.qry_right_overlen = x.qryendpos
            x.start = x.refstartpos - x.qry_left_overlen
            x.end = x.refendpos + x.qry_right_overlen 
            x.ref_left_overlen = x.refstartpos
            x.ref_right_overlen = x.reflen - x.refendpos
            if (x.orientation == '+'):
                x.refstart = x.qrystartpos - x.ref_left_overlen
                x.refend = x.qryendpos + x.ref_right_overlen
            else:
                x.refstart = x.qryendpos - x.ref_right_overlen
                x.refend = x.qrystartpos + x.ref_left_overlen

    num_alms = 0
    for ref in qualify_alms:
        for qry in qualify_alms[ref]:
            if removed[ref, qry] == False:
                num_alms += 1
    current_alms = copy_alms(qualify_alms, removed)
    output_alms(current_alms, output_dir+"/opt_"+str(optmap_i)+"_alms_2_scaled.log")
    print "After scaling, the number of alignments is ", num_alms


    # read qry map
    qry_markers = {}
    with open(myfile+'_r.cmap') as f_q:
        for i in range(11): # 10 lines of header
            header_line = f_q.readline()
        for line in f_q:
            line = line.strip()
            cols = line.split('\t')
            CMapId = int(cols[0])
            ContigLength = float(cols[1])
            NumSites = int(cols[2])       
            SiteID = int(cols[3])
            LabelChannel = cols[4]
            Position = float(cols[5])
            if LabelChannel == "0":
                continue
            if CMapId not in qry_markers:
                qry_markers[CMapId] = []
            Position /= scaling
            qry_markers[CMapId].append(Position)
    for CMapId in qry_markers:
        qry_markers[CMapId].sort()
    f_q.close()

    print '---------------candidate cutting sites-------------------'
    fpair = file(output_dir+"/chimeric_pairs_"+str(optmap_i)+".log", 'w')
    fpair.write("ref_id\tref_pos\tqry_id\tqry_pos\n")
    chimeric_pairs = []
    
    for ref in qualify_alms:
        for qry in qualify_alms[ref]:
            if removed[ref, qry] == True:
                continue
            x = qualify_alms[ref][qry]

            if (x.confidence > min_confidence):
                ref_left_overlen = x.refstartpos
                ref_right_overlen = x.reflen - x.refendpos
                flag_left = False
                flag_right = False
                if (x.qry_left_overlen > minqryoverhang and ref_left_overlen > minrefoverhang and markers_in_qry_left_overhang(qry_markers, x) > 0): 
                    flag_left = True
                    chimeric_pairs.append((x.ref, x.refstartpos, x.qry, x.qrystartpos))
                    print (x.ref, x.refstartpos, x.qry, x.qrystartpos), "is a pair of candidate cutting sites"
                    fpair.write(str(x.ref)+"\t"+str(x.refstartpos)+"\t"+str(x.qry)+"\t"+str(x.qrystartpos)+"\n")
                if (x.qry_right_overlen > minqryoverhang and ref_right_overlen > minrefoverhang and markers_in_qry_right_overhang(qry_markers, x) > 0):
                    flag_right = True
                    chimeric_pairs.append((x.ref, x.refendpos, x.qry, x.qryendpos))
                    print (x.ref, x.refendpos, x.qry, x.qryendpos), "is a pair of candidate cutting sites"
                    fpair.write(str(x.ref)+"\t"+str(x.refendpos)+"\t"+str(x.qry)+"\t"+str(x.qryendpos)+"\n")
                if flag_left == True and flag_right == True:
                    removed[ref, qry] = True
    fpair.close()
    print '---------------END-------------------'
    num_alms = 0
    for ref in qualify_alms:
        for qry in qualify_alms[ref]:
            if removed[ref, qry] == False:
                num_alms += 1
    current_alms = copy_alms(qualify_alms, removed)
    output_alms(current_alms, output_dir+"/opt_"+str(optmap_i)+"_alms_3_removed_both_overhang.log")
    print "After removing alignments with both overhangs, the number of alignments is ", num_alms


    # check overlap between alignments
    for r in qualify_alms:
        for q1 in qualify_alms[r]:
            if removed[r, q1] == True:
                continue
            x = qualify_alms[r][q1]
            for q2 in qualify_alms[r]:
                if removed[r, q2] == True:
                    continue 
                y = qualify_alms[r][q2]
                if q1 >= q2:
                    continue
                if x.refstartpos <= y.refstartpos and y.refstartpos <= x.refendpos:
                    overlap = min(x.refendpos, y.refendpos) - y.refstartpos 
                elif y.refstartpos <= x.refstartpos and x.refstartpos <= y.refendpos:
                    overlap = min(x.refendpos, y.refendpos) - x.refstartpos
                else:
                    overlap = 0
                if overlap >= 20000:
                    if x.confidence < y.confidence:
                        removed[r,q1] = True
                    else:
                        removed[r,q2] = True
    num_alms = 0
    for ref in qualify_alms:
        for qry in qualify_alms[ref]:
            if removed[ref, qry] == False:
                num_alms += 1
    current_alms = copy_alms(qualify_alms, removed)
    output_alms(current_alms, output_dir+"/opt_"+str(optmap_i)+"_alms_4_solved_overlaps.log")
    print "After removing one of two overlap alignments, the number of alignments is ", num_alms
 

    return current_alms, optmap, chimeric_pairs



def core(optmap_list, fasta_file, myfile_list, myfile2_list, output_dir, GLPSOL):
    
    # discard alignments below min_confidence
    #min_confidence = 25
#    header_lines = 10
#    header = []
    # alignment overhangs above this number of bps are considered chimeric
    #minrefoverhang = 100000
    #minqryoverhang = 100000

    min_confidence = params.threshold_1
    false_alm_threshold = params.threshold_2
    num_threads = params.num_threads

    all_alms_list = [] # stores all the Alignments for all groups, all_groups[ref] should contain molecule ref
    all_optmaps = []
    all_chimeric_pairs = []
    for i in range(len(myfile_list)):
        myfile = myfile_list[i]
        myfile2 = myfile2_list[i]
        optmap_file = optmap_list[i]
        alms, optmap,chimeric_pairs = pre_process(i, optmap_file, myfile, myfile2, output_dir, min_confidence)
        all_alms_list.append(alms)
        all_optmaps.append(optmap)
        all_chimeric_pairs.append(chimeric_pairs)

    qualify_alms = {} # only keep one alignment(the one with highest confidence) for each contig in one molecule
    removed = {} # removed[ref,qry] == True means alignment for (ref, qry) is already removed
    qualify_optmap = {}
    
    newref_oldref = {}
    oldref_newref = {}
    new_ref_id = 0
    for i in range(len(myfile_list)):
        alms = all_alms_list[i]
        optmap = all_optmaps[i]
        for ref in optmap:
            new_ref_id += 1
            newref_oldref[new_ref_id] = (i, ref)
            oldref_newref[i, ref] = new_ref_id
            qualify_optmap[new_ref_id] = optmap[ref]
            if ref not in alms:
                continue
            qualify_alms[new_ref_id] = {}
            for qry in alms[ref]:
                x = alms[ref][qry]
                x.ref = new_ref_id
                qualify_alms[new_ref_id][qry] = x
    for newref in qualify_optmap:
        if newref not in qualify_alms:
            qualify_alms[newref] = {}

    fo = open(output_dir+"/newref_oldref.log", 'w')
    for newref in newref_oldref:
        oldref = newref_oldref[newref]
        fo.write(str(newref)+"\t"+str(oldref)+"\n")
    fo.close()


    qualify_chimeric_pairs = []
    chimeric_sites = {}
    for optmap_i in range(len(all_chimeric_pairs)):
        for (ref, refpos, qry, qrypos) in all_chimeric_pairs[optmap_i]:
            newref = oldref_newref[optmap_i, ref]
            qualify_chimeric_pairs.append((newref, refpos, qry, qrypos))
            if newref not in chimeric_sites:
                chimeric_sites[newref] = []
            chimeric_sites[newref].append(refpos)
    for newref in qualify_optmap:
        if newref not in chimeric_sites:
            chimeric_sites[newref] = []
    for newref in chimeric_sites:        
        chimeric_sites[newref].sort()
    

    num_qualify_alms = 0
    for ref in qualify_alms:
        num_qualify_alms += len(qualify_alms[ref])
    print "In total, the number of alignments in qualify_alms is ", num_qualify_alms
    # initialize removed array
    for ref in qualify_alms:
        for qry in qualify_alms[ref]:
            removed[ref,qry] = False


    current_alms = copy_alms(qualify_alms, removed)
    output_alms(current_alms, output_dir+"/alms_0_initial.log")
    print "Initially, the number of alignments is", count_alms(current_alms)    
    alms_0 = copy_alms(qualify_alms, removed)
    aligned_contigs = different_contigs(alms_0, {})
    output_contigs(aligned_contigs, myfile+'_aligned.txt')
    print '---------------END-------------------'


   
    print '---------------removing false positive alignments-------------------'
    current_alms = copy_alms(qualify_alms, removed)
    false_alms(GLPSOL, false_alm_threshold, current_alms, removed, output_dir)   
    current_alms = copy_alms(qualify_alms, removed)
    output_alms(current_alms, output_dir+"/alms_1_removed_false_alms.log")
    print "After removing false positive alignments, the number of alignments is", count_alms(current_alms) 
    print '---------------END-------------------'


    # get repeats
    print '---------------getting repetitive regions------------------'
      # estimate variance
    var = 0.0
    N_var = 0
    for ref in qualify_optmap:
        s, n = get_var(qualify_optmap[ref])
        var += s
        N_var += n
    
    if N_var == 0:
        var = 0
    else:
        var /= N_var
    all_repeats = {}
    c = 0
    for ref in qualify_optmap:
        c += 1
        repeats = get_repeats(qualify_optmap[ref], var, num_threads)
        all_repeats[ref] = repeats
        if ref not in qualify_alms:
            continue
    print '---------------END-------------------'




       
    print '---------------removing contained contigs locally-------------------'
    for ref in qualify_alms:
        for q1 in qualify_alms[ref]:
            x = qualify_alms[ref][q1]
            for q2 in qualify_alms[ref]:
                if q2 <= q1:
                    continue
                y = qualify_alms[ref][q2]
                if (x.start >= y.start) and (x.end <= y.end):
                    removed[ref, q1] = True
                    print [ref, q1], "alignment is removed becasue it's contained in alignment", [ref, q2]
                elif (y.start >= x.start) and (y.end <= x.end):
                    removed[ref, q2] = True
                    print [ref, q2], "alignment is removed becasue it's contained in alignment", [ref, q1]
    current_alms = copy_alms(qualify_alms, removed)
    output_alms(current_alms, output_dir+"/alms_2_removed_contained_locally.log")
    print "After removing contained alignments locally, the number of alignments is", count_alms(current_alms)
    print '---------------END-------------------'



    #build new mst
    print '---------------building new mst-------------------'
    fo = file(output_dir+"/ugraph.log", 'w')
    current_alms = copy_alms(qualify_alms, removed) 
    forest, vertex_orientations = get_mst(current_alms, fo)
    fo.close()
    output_forest(forest, vertex_orientations, output_dir+"/forest.log")
    print '---------------END-------------------'


    # build order graph
    current_alms = copy_alms(qualify_alms, removed)
    order_graphs, contigs_oris, raw_gaps_list = merge_order_graphs(current_alms, all_repeats, chimeric_sites, forest, vertex_orientations)
    output_DAGs(order_graphs, output_dir+"/ordergraph.log")


    fo = file(output_dir+"/contig_oris.txt", 'w')
    for qry in contigs_oris:
        fo.write(str(qry)+"\t"+str(contigs_oris[qry])+"\n")
    fo.close()

    # MinULP solving
    order_graphs_fullorder = [] 
    for i in range(len(order_graphs)):
        og = order_graphs[i]
        og = minFAS(og)
        og_unsorted = []
        for v in og:
            edges = list(og[v])
            og_unsorted.append((v, edges))
        if not check_full_order(og_unsorted):
            candiEdges = get_candiedges(og)
            if len(candiEdges) <= 20:
                og_new_set = full_order.exhaust(og)
            else:
                og_new_set = full_order.apprx(og,GLPSOL, output_dir, i)
            for og_new in og_new_set:
                order_graphs_fullorder.append(og_new)
        else:
            order_graphs_fullorder.append(og)

    orders = []
    for og in order_graphs_fullorder:
        og_unsorted = []
        for v in og:
            edges = list(og[v])
            og_unsorted.append((v, edges))
        og_sorted = topolgical_sort(og_unsorted) 
        order = []
        for item in og_sorted:
            order.append(item[0])
        order = list(reversed(order))
        orders.append(order)

    fo = file(output_dir+"/orders.txt", 'w')
    for order in orders:
        for item in order:
            fo.write(str(item)+"\t")
        fo.write("\n")
    fo.close()

    # gap estimation
    gaps_list = []
    Edges = set([])
    for og in order_graphs_fullorder:
        for v in og:
            for u in og[v]:
                Edges.add((v,u))
    for v, u, gap in raw_gaps_list:
        if (v, u) in Edges:
            gaps_list.append((v, u, gap))

    contigs_lens = {}
    num = 1
    for line in open(fasta_file):
        line = line.strip()
        if line[0] == '>':
            current_id = num
            num += 1
            contigs_lens[current_id] = 0
        else:
            contigs_lens[current_id] += len(line)

    final_gaps = gaps(orders, gaps_list, contigs_lens)
    
    fo = file(output_dir+"/gaps.txt", 'w')
    for qry1, qry2 in final_gaps:
        fo.write(str(qry1)+"\t"+str(qry2)+"\t"+str(final_gaps[qry1, qry2])+"\n")
    fo.close()



