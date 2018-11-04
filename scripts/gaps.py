import numpy as np
from scipy.optimize import minimize
from scipy.special import loggamma
import math
from params import gap_naive

def func_log_1(params, obs):
    k = sum(params)
    return (k/2-1)*np.log(obs)-obs/2-(k/2)*np.log(2)-loggamma(k/2).real
    
def func_logl(params, samples):
    logl = 0
    for obs, begin, end in samples:
        logl += func_log_1(params[begin:end+1], obs)
    return (-logl)

def gaps(orders, gaps_list, qry_lens):
    # naive method to get initial values
    initial_gaps = {}
    for order in orders:
        for i in range(0, len(order)-1):
            qry1 = order[i]
            qry2 = order[i+1]
            gap_sum = 0.0
            gap_num = 0
            for (u, v, gap) in gaps_list:
                if (qry1, qry2) == (u, v):
                    gap_sum += gap
                    gap_num += 1
            gap_ave = gap_sum / gap_num
            initial_gaps[qry1,qry2] = gap_ave
    if gap_naive == True:
        return initial_gaps

    # get qry_index 
    qry_index = {}
    for i in range(len(orders)):
        order = orders[i]
        for j in range(len(order)):
            qry = order[j]
            qry_index[qry] = (i,j)
    # clustring
    gaps_cluster = []
    for order in orders:
        gaps_cluster.append([i for i in range(len(order)-1)])

          # deal with negative gap
    for qry1, qry2 in initial_gaps:
        (i_1, j_1) = qry_index[qry1]
        if initial_gaps[qry1, qry2] < 0:
            gaps_cluster[i_1][j_1] = -1 
    gaps_list_new = []
    for (u, v, gap) in gaps_list:
        (u_i, u_j) = qry_index[u]
        (v_i, v_j) = qry_index[v]
        if u_i != v_i or u_j >= v_j:
            print("ERROR!!!")  
            exit()
        # check if this contig contains negative gap
        neggap = False
        for k in range(u_j, v_j):
            if gaps_cluster[u_i][k] == -1:
                neggap = True
                break
        if neggap == True:
            continue
        gaps_list_new.append((u, v, gap))
        # clustering
        new_cluster_id = gaps_cluster[u_i][u_j]
        last_cluster_id = gaps_cluster[u_i][v_j-1]
        for k in range(u_j, v_j):
            gaps_cluster[u_i][k] = new_cluster_id
        for k in range(v_j, len(gaps_cluster[u_i])):
            if gaps_cluster[u_i][k] != last_cluster_id:
                break
            gaps_cluster[u_i][k] = new_cluster_id
 
    # grouping samples to each cluster
    samples = {}
    for (u, v, gap) in gaps_list_new:
        (u_i, u_j) = qry_index[u]
        cluster_id = gaps_cluster[u_i][u_j]
        if (u_i, cluster_id) not in samples:
            samples[u_i, cluster_id] = []
        samples[u_i, cluster_id].append((u, v, gap))


    
    # remove contig length in gap distances
    samples_trsf = {}
    for (s_id, cluster_id) in samples:
        slist = []
        for i in range(len(samples[s_id, cluster_id])):
            (u, v, gap) = samples[s_id, cluster_id][i]
            (u_i, u_j) = qry_index[u]
            (v_i, v_j) = qry_index[v]
            for k in range(u_j+1, v_j):
                qry = orders[u_i][k]
                qry_len = qry_lens[qry]
                gap -= qry_len
            slist.append((gap, u_j-cluster_id, v_j-cluster_id-1))    
        samples_trsf[s_id, cluster_id] = slist

    # estimation
    final_gaps = {}
    for (s_id, cluster_id) in samples_trsf:
        sub_sample = samples_trsf[s_id, cluster_id]
        x0 = []
        for k in range(len(gaps_cluster[s_id])):
            if cluster_id == gaps_cluster[s_id][k]:
                qry1 = orders[s_id][k]
                qry2 = orders[s_id][k+1]
                x0.append(initial_gaps[qry1,qry2])
        flagneg = False 
        for obs, begin, end in sub_sample:
            if obs <= 0:
                flagneg = True
        if flagneg == True:
            for k in range(len(x0)):
                qry1 = orders[s_id][cluster_id + k]
                qry2 = orders[s_id][cluster_id + k + 1]
                final_gaps[qry1, qry2] = initial_gaps[qry1, qry2]              
            continue            
        res = minimize(func_logl, x0, args=(sub_sample,), method='BFGS')
 
        params = res.x
        for k in range(len(params)):
            qry1 = orders[s_id][cluster_id + k]
            qry2 = orders[s_id][cluster_id + k + 1]
            final_gaps[qry1, qry2] = params[k]              

    # deal with negative gap 
    for qry1, qry2 in initial_gaps:
        if initial_gaps[qry1, qry2] < 0:
            final_gaps[qry1,qry2] = initial_gaps[qry1, qry2]            
    return final_gaps

