from graph.undirected_graph import UndirectedGraph
from util.permutation import selection_sort
from sage.graphs.graph import Graph
from sage.graphs.graph_generators import GraphGenerators

def nauty_canonicalize(g):
    n = len(g)
    edges = g.edges()
    G, sigma = Graph([list(range(n)), edges]).canonical_label(certificate=True)
    new_edges = list(G.edges(labels=False))
    edge_permutation = [tuple(sorted([sigma[edge[0]],sigma[edge[1]]])) for edge in edges]
    index_permutation = [new_edges.index(e) for e in edge_permutation]
    return UndirectedGraph(n, list(new_edges)), selection_sort(index_permutation)

def nauty_has_odd_automorphism(g):
    n = len(g)
    edges = g.edges()
    G = Graph([list(range(n)), edges])
    for sigma in G.automorphism_group().gens():
        edge_permutation = [tuple(sorted([sigma(edge[0]),sigma(edge[1])])) for edge in edges]
        index_permutation = [edges.index(e) for e in edge_permutation]
        if selection_sort(index_permutation) == -1:
            return True
    return False

def nauty_generate(num_vertices, num_edges, connected=None, biconnected=None, min_degree=0):
    args = "{} {}:{}".format(num_vertices, num_edges, num_edges)
    if connected:
        args += " -c"
    if biconnected:
        args += " -C"
    if min_degree != 0:
        args += " -d{}".format(min_degree)
    for G in GraphGenerators().nauty_geng(args):
        G = G.canonical_label()
        yield UndirectedGraph(num_vertices, list(G.edges(labels=False)))
