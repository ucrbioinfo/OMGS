def minFAS(G): # greedy algorithm
    g = G.copy()
    remove_edge_set = set([])
    while True:
        SCC = split2SCC(g)
        g = simplify(SCC, g)        
        if g == {}:
            break
        e_set = remove_edges(g)
        for e in e_set:
            remove_edge_set.add(e)

    # remove edges from G
    new_G = {}
    for v in G:
        new_G[v] = {}
        for u in G[v]:
            if (v, u) not in remove_edge_set:
                new_G[v][u] = G[v][u]
    return new_G

def simplify(SCC, g):
    mapv2u = {}
    for u in SCC:
        scc = SCC[u]
        for v in scc:
            mapv2u[v] = u
    new_g = {}
    for v in g:
        u =  mapv2u[v] 
        if len(SCC[u]) == 1:
            continue 
        new_g[v] = {}
        for k in g[v]:
            if k not in SCC[u]:
                continue
            new_g[v][k] = g[v][k]
    return new_g
    

def remove_edges(g):
    indegree = {}
    outdegree = {}
    score = {}
    for v in g:
        indegree[v] = 0
        outdegree[v] = 0
    for v in g:
        for u in g[v]:
            indegree[u] += g[v][u]
            outdegree[v] += g[v][u]
    for v in g:
        score[v] = max(indegree[v]/outdegree[v], outdegree[v]/indegree[v])
    maxs = -1
    for v in g:
        if score[v] > maxs:
            maxs = score[v]
            maxv = v
    remove_edge_set = set([])
    if indegree[maxv] > outdegree[maxv]:
        # remove out edges
        for u in g[maxv]:
            remove_edge_set.add((maxv,u))
        g[maxv] = {} 
    else:
        # remove in edges
        for u in g:
            for v in g[u]:
                if v != maxv:
                    continue
                remove_edge_set.add((u,v))
        for (u,v) in remove_edge_set:
            del g[u][v]      
    return remove_edge_set

def split2SCC(g): # find strong connected component
    f = {}
    DFS(g, f) 
    g_t = transpose(g)
    f_sorted = sorted(f.items(), key=lambda x: x[1], reverse=True)
    SCC = DFS2(g_t, f_sorted)

    return SCC

def transpose(g):
    g_t = {}
    for v in g:
        g_t[v] = {}
    for v in g:
        for u in g[v]:
            g_t[u][v] = g[v][u]
    return g_t 

def DFS(g, f):
    color = {}
    for u in g:
        color[u] = 0 # white
    time = [0]
    for u in g:
        if color[u] == 0:
            DFS_visit(g, u ,time, color, f)

def DFS_visit(g, u, time, color, f):
    time[0] += 1
    color[u] = 1 # gray
    for v in g[u]:
        if color[v] == 0:
            DFS_visit(g, v, time, color, f)
    color[u] = 2
    time[0] += 1
    f[u] = time[0]


def DFS2(g, f_sorted):
    color = {}
    for u in g:
        color[u] = 0 # white
    time = [0]
    SCC = {}
    for (u,value) in f_sorted:
        if color[u] == 0:
            SCC[u] = set([])
            DFS_visit2(g, u ,time, color, SCC[u])
    return SCC

def DFS_visit2(g, u, time, color, scc):
    scc.add(u)
    time[0] += 1
    color[u] = 1 # gray
    for v in g[u]:
        if color[v] == 0:
            DFS_visit2(g, v, time, color, scc)
    color[u] = 2
    time[0] += 1




G = {}
G['a'] = {}
G['b'] = {}
G['c'] = {}
G['d'] = {}
G['e'] = {}
G['f'] = {}
G['g'] = {}
G['h'] = {}
G['a']['b']=1
G['b']['c']=2
G['b']['e']=3
G['b']['f']=4
G['c']['d']=5
G['c']['g']=6
G['d']['c']=7
G['d']['h']=8
G['e']['a']=9
G['e']['f']=10
G['f']['g']=11
G['g']['f']=12
G['g']['h']=13
#G['h']['h']=14

