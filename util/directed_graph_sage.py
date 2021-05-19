from graph.directed_graph import DirectedGraph
from util.permutation import selection_sort
from itertools import combinations
import sage.all # make SageMath work when called from Python
from sage.graphs.digraph import DiGraph
from sage.graphs.graph import Graph
from sage.graphs.graph_generators import GraphGenerators

# TODO: use directg?
# TODO: cache

def directed_graph_canonicalize(g):
    n = len(g)
    edges = g.edges()
    G, sigma = DiGraph([list(range(n)), edges]).canonical_label(certificate=True)
    new_edges = list(G.edges(labels=False))
    edge_permutation = [tuple([sigma[edge[0]],sigma[edge[1]]]) for edge in edges]
    index_permutation = [new_edges.index(e) for e in edge_permutation]
    undo_canonicalize = [0]*n
    for k, v in sigma.items():
        undo_canonicalize[v] = k
    return DirectedGraph(n, list(new_edges)), undo_canonicalize, selection_sort(index_permutation)

def directed_graph_has_odd_automorphism(g):
    n = len(g)
    edges = g.edges()
    G = DiGraph([list(range(n)), edges])
    for sigma in G.automorphism_group().gens(): # NOTE: it suffices to check generators
        edge_permutation = [tuple([sigma(edge[0]),sigma(edge[1])]) for edge in edges]
        index_permutation = [edges.index(e) for e in edge_permutation]
        if selection_sort(index_permutation) == -1:
            return True
    return False

def directed_graph_generate(num_vertices, num_edges, connected=None, biconnected=None, min_degree=0, loops=True):
    extra_args = ""
    if connected:
        extra_args += " -c"
    if biconnected:
        extra_args += " -C"
    if min_degree != 0:
        extra_args += " -d{}".format(min_degree)
    try:
        if loops:
            max_loop_order = num_edges // 2 # NOTE: can have at most this many loops, while still attaining num_edges
        else:
            max_loop_order = 0
        for loop_order in range(max_loop_order + 1):
            nauty_args = "{} {}:{} {}".format(num_vertices, num_edges - loop_order, num_edges - loop_order, extra_args)
            for G in GraphGenerators().nauty_geng(nauty_args):
                results = []
                for loops in combinations(G.edges(), loop_order): # NOTE: no need for triple edges, they always yield zero graphs
                    H = Graph([G.vertices(), G.edges()], multiedges=True)
                    H.add_edges(loops)
                    for h in H.orientations():
                        if h.has_multiple_edges(): # NOTE: graphs with directed multiple edges are zero
                            continue
                        h = h.canonical_label()
                        g = DirectedGraph(num_vertices, list(h.edges(labels=False)))
                        g.canonicalize_edges()
                        if not g in results:
                            yield g
                            results.append(g)
    except ValueError: # impossible values also errors
        pass
