from graph.undirected_graph import UndirectedGraph
from util.permutation import selection_sort
import sage.all # make SageMath work when called from Python
from sage.graphs.graph import Graph
import subprocess
import os

NAUTY_PREFIX = '' # e.g. '/home/rburing/src/nauty27r1/'

def nauty_generate_undirected(num_vertices, num_edges, connected=None, biconnected=None, min_degree=0):
    args = [str(num_vertices), "{}:{}".format(num_edges, num_edges)]
    if connected:
        args.append("-c")
    if biconnected:
        args.append("-C")
    if min_degree != 0:
        args.append("-d{}".format(min_degree))
    FNULL = open(os.devnull, 'w')
    geng = subprocess.Popen((NAUTY_PREFIX + 'geng', *args), stdout=subprocess.PIPE, stderr=FNULL)
    showg = subprocess.Popen((NAUTY_PREFIX + 'showg', '-e', '-l0'), stdin=geng.stdout, stderr=FNULL, stdout=subprocess.PIPE)
    line_count = -1
    for line in showg.stdout:
        if line_count % 4 == 2:
            graph_encoding = line.decode('ascii').rstrip()
            edges = [tuple(map(int,e.split(' '))) for e in graph_encoding.split('  ')]
            yield Graph([list(range(num_vertices)), edges])
        line_count += 1

def undirected_graph_canonicalize(g):
    n = len(g)
    edges = g.edges()
    G, sigma = Graph([list(range(n)), edges]).canonical_label(certificate=True)
    new_edges = list(G.edges(labels=False))
    edge_permutation = [tuple(sorted([sigma[edge[0]],sigma[edge[1]]])) for edge in edges]
    index_permutation = [new_edges.index(e) for e in edge_permutation]
    undo_canonicalize = [0]*n
    for k, v in sigma.items():
        undo_canonicalize[v] = k
    return UndirectedGraph(n, list(new_edges)), undo_canonicalize, selection_sort(index_permutation)

def undirected_graph_has_odd_automorphism(g):
    n = len(g)
    edges = g.edges()
    G = Graph([list(range(n)), edges])
    for sigma in G.automorphism_group().gens(): # NOTE: it suffices to check generators
        edge_permutation = [tuple(sorted([sigma(edge[0]),sigma(edge[1])])) for edge in edges]
        index_permutation = [edges.index(e) for e in edge_permutation]
        if selection_sort(index_permutation) == -1:
            return True
    return False

def undirected_graph_generate(num_vertices, num_edges, connected=None, biconnected=None, min_degree=0, has_odd_automorphism=None):
    for G in nauty_generate_undirected(num_vertices, num_edges, connected=connected, biconnected=biconnected, min_degree=min_degree):
        G = G.canonical_label()
        g = UndirectedGraph(num_vertices, list(G.edges(labels=False)))
        if has_odd_automorphism is None or undirected_graph_has_odd_automorphism(g) == has_odd_automorphism:
            yield g
