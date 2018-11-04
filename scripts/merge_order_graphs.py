#!/usr/bin/python
import csv
import sys

from sets import Set

from collections import defaultdict

import collections
import math
from math import exp

def DFS_DAG(ug, g, visited, start, subgraph):
    visited[start] = True
    subgraph[start] = g[start]
    preds = ug[start]
    for v in preds:
        if visited[v] == True:
            continue
        DFS_DAG(ug, g, visited, v, subgraph)

def splite_graph(g):
    ug = {}
    for v in g:
        ug[v] = set([])
    for v in g:
        for u in g[v]:
            ug[v].add(u)
            ug[u].add(v)
    visited = {}
    for v in g:
        visited[v] = False
    new_g_set = []
    while True:
        found = False
        for v in g:
            if visited[v] == False:
                found = True
                start = v
                break
        if found == True:
            subgraph = {}
            DFS_DAG(ug, g, visited, start, subgraph)
            new_g_set.append(subgraph)
        else:
            break
    return new_g_set

def splite_graphs(graph_set):
    new_graph_list = []
    for root in graph_set:
        g = graph_set[root]
        new_g_set = splite_graph(g)
        for new_g in new_g_set:
            new_graph_list.append(new_g)
    return new_graph_list

def DFS(g, v, u, visited):#search for path from u to v
    if u == v:
        return True
    visited[u] = True
    for w in g[u]:
        if visited[w] == True:
            continue
        if DFS(g, v, w, visited) == True:
            return True
    return False       

def merge_DAG(g1,oris_1,  g2, oris_2):
    for v in g2:
        if v not in g1:
            g1[v] = {}
            oris_1[v] = oris_2[v]
    for v in g2:
        next_nodes = g2[v]
        for u in next_nodes:
            if u in g1[v]:
                g1[v][u] += g2[v][u]
                continue
            g1[v][u] = g2[v][u]



def contain(query, start, end):
    if start > end:
        return False
    if query >= start and query <= end:
        return True
    else:
        return False

def get_optmap_quality(begin, end, repeats_ref, chimeric_sites_ref):
    for site in chimeric_sites_ref:
        if abs(begin-site) < 5000 or abs(end-site) < 5000:
            return -1
    num_chimeric = 0
    for site in chimeric_sites_ref:
        if contain(site, begin-10000, end+10000):
            num_chimeric += 1
    if num_chimeric > 0:
        return -1


    if begin >= end:
        print "Note: begin is larger than end!!!"
        repeats_len = 0
        return exp(-1.0*(repeats_len/100000+1))        

    repeats_len = 0.0
    for (r_start, r_end) in repeats_ref:
        if contain(begin, r_start, r_end) and contain(end, r_start, r_end):
            repeats_len += (end - begin)
        elif contain(begin, r_start, r_end) and not contain(end, r_start, r_end):
            repeats_len += (r_end - begin)
        elif not contain(begin, r_start, r_end) and contain(end, r_start, r_end):
            repeats_len += (end - r_start)
        elif begin < r_start and r_end < end:
            repeats_len += r_end - r_start
    
    optmap_quality = exp(-1.0*(repeats_len/100000+1))
    return optmap_quality

def build_DAG2(alms, repeats, chimeric_sites, ref, ori, raw_gaps_list):
    DAG = {}
    for qry in alms[ref]:
        DAG[qry] = {}
    vertex_oris = {}
    sorted_qrys = []
    for qry in alms[ref]:
        DAG[qry] = {}
        sorted_qrys.append((qry, alms[ref][qry].refstartpos))
    if ori == 1:
        sorted_qrys.sort(key=lambda x: x[1])
    elif ori == -1:
        sorted_qrys.sort(key=lambda x: x[1], reverse=True)
    for qry in alms[ref]:
        if ori == 1:
            vertex_oris[qry] = alms[ref][qry].orientation
        elif ori == -1:
            if alms[ref][qry].orientation == "+":
                vertex_oris[qry] = '-'
            elif alms[ref][qry].orientation == "-":
                vertex_oris[qry] = '+'

    for i in range(0, len(sorted_qrys)-1):
        (qry1, start1) = sorted_qrys[i]
        (qry2, start2) = sorted_qrys[i+1]
        if qry1 == 21 and qry2 == 758:
            print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
            print ori
            print ref
            print chimeric_sites[ref]
            print alms[ref][qry1].refendpos, alms[ref][qry2].refstartpos
        if ori == 1:
            optmap_quality = get_optmap_quality(alms[ref][qry1].refendpos, alms[ref][qry2].refstartpos, repeats[ref], chimeric_sites[ref])
            gap = alms[ref][qry2].start - alms[ref][qry1].end
        elif ori == -1:
            optmap_quality = get_optmap_quality(alms[ref][qry2].refendpos, alms[ref][qry1].refstartpos, repeats[ref], chimeric_sites[ref])
            gap = alms[ref][qry1].start - alms[ref][qry2].end
        if qry1 == 21 and qry2 == 758:
            print optmap_quality
        if optmap_quality > 0:
            DAG[qry1][qry2] = optmap_quality * (alms[ref][qry1].confidence + alms[ref][qry2].confidence)
            raw_gaps_list.append((qry1, qry2, gap))
    return DAG, vertex_oris

def tree_traversal_merge_DAGs(alms, repeats, chimeric_sites,onetree, onetree_v_oris, root, DAG, vertex_oris,  raw_gaps_list):
    this_ori = onetree_v_oris[root]
    DAG_one, vertex_oris_one = build_DAG2(alms, repeats, chimeric_sites, root, this_ori, raw_gaps_list)
    merge_DAG(DAG, vertex_oris, DAG_one, vertex_oris_one)
    #search children
    children = onetree[root] 
    for c in children:
        tree_traversal_merge_DAGs(alms, repeats, chimeric_sites, onetree, onetree_v_oris, c, DAG, vertex_oris, raw_gaps_list) 

def merge_order_graphs(alms, repeats, chimeric_sites, forest, vertex_orientations):
    #merge order graphs
    
    raw_gaps_list = []

    DAG_set = {}
    vertex_oris_set= {}
    for root in forest:
        DAG = {}
        vertex_oris = {}
        onetree = forest[root]
        onetree_v_oris = vertex_orientations[root]
        tree_traversal_merge_DAGs(alms, repeats, chimeric_sites,  onetree, onetree_v_oris, root, DAG, vertex_oris, raw_gaps_list)
        DAG_set[root] = DAG
        vertex_oris_set[root] = vertex_oris
    
    vertex_oris_all = {}
    for root in DAG_set:
        for v in DAG_set[root]:
            if v not in vertex_oris_all:
                vertex_oris_all[v] = vertex_oris_set[root][v]


    DAG_list= splite_graphs(DAG_set)
    

    print "In total, the number of connected DAG is", len(DAG_list)

    #draw
    f = file("order_graph.dot",'w')
    f.write("digraph {\n")
    for d in DAG_list:
        for v in d:
            for u in d[v]:
                f.write(str(v)+" -> "+str(u)+"[label=\""+str(d[v][u])+"\"];\n")
    f.write("}\n")
    return DAG_list, vertex_oris_all, raw_gaps_list

