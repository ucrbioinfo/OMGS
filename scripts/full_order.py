import copy

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



def check_full_order(graph_unsorted):
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
        candidate_set = []
        for node, edges in list(graph_unsorted.items()):
            for edge in edges:
                if edge in graph_unsorted:
                    break
            else:
                acyclic = True
                candidate_set.append(node)

        if not acyclic:
            return False
        if len(candidate_set) > 1:
            return False

        node = candidate_set[0]
        del graph_unsorted[node]

    return True

def search(g, Edges, opt_solution, opt_value, solution, m, n):
    if m == n:
        #check if solution is feasible
        value = 0
        g_cp = copy.deepcopy(g) 
        for i in range(0, len(solution)):
            if solution[i] == False:
                continue
            (v, u) = Edges[i]
            value += g[v][u]
            # remove (v, u) from g_cp
            del g_cp[v][u]
        feasible = True
        g_cp_set = splite_graph(g_cp)
        for g_sub in g_cp_set:
            g_sub_unsorted = []
            for v in g_sub:
                edges = list(g_sub[v])
                g_sub_unsorted.append((v, edges))
            if check_full_order(g_sub_unsorted) == False:
                feasible = False
                break 

        if feasible == False:
            return

        #compare with current best solution
        if value < opt_value[0]:
            for i in range(0, len(solution)):
                opt_solution[i] = solution[i]
            opt_value[0] = value
            
    else:
        solution[m] = False
        search(g, Edges, opt_solution, opt_value, solution, m+1, n)
        solution[m] = True
        search(g, Edges, opt_solution, opt_value, solution, m+1, n)

def get_candiedges(g):
    g_trans = {}
    for v in g:
        g_trans[v] = {}
    for v in g:
        for u in g[v]:
            g_trans[u][v] = g[v][u]
    conjunction = set([])
    for v in g:
        if len(g[v]) >1 or len(g_trans[v]) >1:
            conjunction.add(v)
    
    candiEdges = []    
    for v in g:
        for u in g[v]:
            if v in conjunction or u in conjunction:
                candiEdges.append((v,u))
    return candiEdges

def exhaust(g):#exhaust Agortihm

    candiEdges = get_candiedges(g)

    solution = []
    opt_solution = []
    for i in range(0, len(candiEdges)):
        solution.append(False)
        opt_solution.append(False)
    opt_value = [0]
    for v in g:
        for u in g[v]:
            opt_value[0] += g[v][u]
    opt_value[0] += 1.0

    search(g, candiEdges, opt_solution, opt_value, solution, 0, len(candiEdges))

    g_cp = copy.deepcopy(g) 
    for i in range(0, len(opt_solution)):
        if opt_solution[i] == False:
            continue
        (v, u) = candiEdges[i]
        # remove (v, u) from g_cp
        del g_cp[v][u]
    g_cp_set = splite_graph(g_cp)

    
    return g_cp_set



def apprx(g,GLPSOL, output_dir, i):#approximation algorithm
    vertices = set([])
    for v in g:
        vertices.add(v)
    reachable = DFS_reach(g)
    edges = {}
    for u in reachable:
        for v in reachable[u]:
            if v in g[u]:
                edges[(min(u,v),max(u,v))] = g[u][v]
            else:
                edges[(min(u,v),max(u,v))] = 0
    ECP(GLPSOL, output_dir, i, vertices, edges)    

def DFS_reach(g):
    reachable = {}
    for u in g:
        reachable[u] = {}
        color = {}
        for v in g:
            color[v] = 0 # white
        DFS_visit(g, u, color, reachable[u])
    return reachable

def DFS_visit(g, u, color, r):
    r.add(u)
    color[u] = 1 # gray
    for v in g[u]:
        if color[v] == 0:
            DFS_visit(g, v, color, r)
    color[u] = 2



