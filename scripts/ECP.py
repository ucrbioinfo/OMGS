import math
def ECP(GLPSOL, output_dir, i, vertices, edges):
    W = max(edges.values())
    new_edges = {}
    for u in vertices:
        for v in vertices:
            if u >= v:
                continue
            if (u, v) in edges:
                new_edges[(u,v)] = edges[(u,v)]
            else:
                new_edges[(u,v)] = -W * ((n * math.log2(n))**2)
    #clusters = corre_clustering(vertices, new_edges)
    clusters = corre_clustering(GLPSOL, output_dir, i, vertices, new_edges)
    return clusters 


