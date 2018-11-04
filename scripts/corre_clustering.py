import os
def create_LP_file(filename, vertices, edges):
    fo = file(filename, 'w')

    # write objective function
    fo.write("Minimize\n")
    fo.write("value: ")
    num = 0
    edge_to_id = {}
    for e in edges:
        num += 1
        edge_to_id[e] = num 
        if edges[e] < 0:
            m = abs(edges[e])
        else:
            m = 0
        if edges[e] > 0:
            p = abs(edges[e])
        else:
            p = 0
        mp = p - m
        if mp >= 0:
            sign = " + "
        else:
            sign = " - "
        fo.write(sign  + str(abs(mp)) + " x" +str(num)   + "\n")
    fo.write("\n")

    # write constrains
    fo.write("Subject To\n")
    vertex_list = list(vertices)
    vertex_list.sort()
    num_constrains = 0
    for i in range(len(vertex_list)):
        u = vertex_list[i]
        for j in range(i+1, len(vertex_list)):
            v = vertex_list[j]
            for k in range(j+1, len(vertex_list)):
                w = vertex_list[k]
                id1 = edge_to_id[(u,v)]
                id2 = edge_to_id[(u,w)]
                id3 = edge_to_id[(v,w)]
                fo.write("x"+str(id1)+" + x"+str(id2)+" - x"+str(id3) +" >= "+str(0)+"\n")
                fo.write("x"+str(id1)+" + x"+str(id3)+" - x"+str(id2) +" >= "+str(0)+"\n")
                fo.write("x"+str(id2)+" + x"+str(id3)+" - x"+str(id1) +" >= "+str(0)+"\n")
                num_constrains += 3
    # write bounds
    fo.write("\n")
    fo.write("Bounds\n")
    edge_list = []
    for e in edges:
        edge_list.append(e)
        id1 = edge_to_id[e]
        fo.write(str(0) +" <= x"+str(id1) +" <= "+str(1)+"\n")
    fo.write("\n")

    fo.write("End\n")
    fo.close()

    return edge_list, num_constrains

def read_results(filename, edge_list, num_constrains):
    LP_values = {}
    with open(filename) as fi:
        fi.readline()
        line = fi.readline()
        line = line.strip()
        cols = line.split()
        opt_value = float(cols[2])
        for i in range(num_constrains): 
            fi.readline()
        for i in range(len(edge_list)):
            e = edge_list[i] 
            line = fi.readline()
            line = line.strip()
            cols = line.split()
            value = float(cols[1])
            LP_values[e] = value
    return LP_values

def update_G(G, B):
    for v in B:
        neighbors = G[v].keys()
        for u in neighbors:
            del G[u][v]
        del G[v]
            
def cut(S, G, edges):
    c = 0
    for v in S:
        for u in G[v]:
            if u in S:
                continue
            e = (min(u,v), max(u,v))
            if edges[e] > 0:
                p = abs(edges[e])
            else:
                p = 0
            c += p
    return c

def vol(S, G, edges):
    c = 0
    for v in S:
        for u in G[v]:
            if u not in S:
                continue
            if u <= v:
                continue
            e = (min(u,v), max(u,v))
            if edges[e] > 0:
                p = abs(edges[e])
            else:
                p = 0
            c += p * G[v][u]
    return c

def cal_B(G, u, r):
    B = set([])
    for v in G[u]:
        if G[u][v] <= r:
            B.add(v)
    return B


def grour(r, u, G, B):
    mind = -1
    
    for v in G:
        if v in B:
            continue
        if G[u][v] <= r:
            continue
        if G[u][v] > mind:
            mind = G[u][v]
    if mind < 0:
        sys.exit("mind cannot be smaller than 0!!")
    newr = mind 
  
    return newr

def round(vertices, edges, LP_values):
    clusters = []
    G = {}
    for v in vertices:
        G[v] = {}
    for (u,v) in LP_values:
        G[u][v] = LP_values[(u,v)]
        G[v][u] = LP_values[(u,v)]

    while(True):
        u = G.keys()[0]
        r = 0
        B = set([])
        while(True):
            r += growr(r, G, B) 
            B = calc_B(G, u, r) 
            if cut(B, G, edges) <= c * math.log(n+1) * vol(B, G, edges):
                break
        clusters.append(B)
        update_G(G, B)
        if G == {}:
            break

    return clusters 

def corre_clustering(GLPSOL, output_dir, i, vertices, edges):#LP based Agortihm
    lp_file = output_dir + "/lp_"+str(i)+".lp"
    lp_result = output_dir + "/lp_"+str(i)+".sol"
    edge_list, num_constrains = create_LP_file(lp_file, vertices, edges)
    command = GLPSOL + " --lp " + lp_file + " -w " + lp_result
    os.system(command)
    LP_values = read_results(lp_result, edge_list, num_constrains)
    print(LP_values)
    clusters = round(vertices, edges, LP_values)
    return clusters

GLPSOL = "glpsol"
vertices = set([3, 4, 5, 6])
edges = {}
edges[(3,4)] = -4
edges[(3,5)] = 3
edges[(3,6)] = -7
edges[(4,5)] = -2
edges[(4,6)] = 8
edges[(5,6)] = 2
corre_clustering(GLPSOL, ".", 0, vertices, edges)





