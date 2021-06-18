from graph.directed_graph import DirectedGraph
from util.permutation import selection_sort
from itertools import combinations
import sage.all # make SageMath work when called from Python
from sage.graphs.digraph import DiGraph
import subprocess
import os

NAUTY_PREFIX = '' # e.g. '/home/rburing/src/nauty27r1/'

def nauty_generate_directed(num_vertices, num_undirected_edges, num_directed_edges, connected=None, biconnected=None, min_degree=0, loops=True):
    geng_args = [str(num_vertices), "{}:{}".format(num_undirected_edges, num_undirected_edges)]
    if connected:
        geng_args.append("-c")
    if biconnected:
        geng_args.append("-C")
    if min_degree != 0:
        geng_args.append("-d{}".format(min_degree))
    directg_args = ['-T', '-e{}'.format(num_directed_edges)]
    if not loops:
        directg_args.append('-o')
    FNULL = open(os.devnull, 'w')
    geng = subprocess.Popen((NAUTY_PREFIX + 'geng', *geng_args), stdout=subprocess.PIPE, stderr=FNULL)
    directg = subprocess.Popen((NAUTY_PREFIX + 'directg', *directg_args), stdin=geng.stdout, stdout=subprocess.PIPE, stderr=FNULL)
    for line in directg.stdout:
        graph_encoding = line.decode('ascii').rstrip()
        numbers = [int(v) for v in graph_encoding.split(' ')]
        numbers = numbers[2:]
        edges = [(numbers[2*i], numbers[2*i+1]) for i in range(num_directed_edges)]
        yield DiGraph([list(range(num_vertices)), edges])

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

def directed_graph_generate(num_vertices, num_edges, connected=None, biconnected=None, min_degree=0, loops=True, has_odd_automorphism=None):
    if loops:
        max_loop_order = num_edges // 2 # NOTE: can have at most this many loops, while still attaining num_edges
    else:
        max_loop_order = 0
    for loop_order in range(max_loop_order + 1):
        for h in nauty_generate_directed(num_vertices, num_edges - loop_order, num_edges, connected=connected, biconnected=biconnected, min_degree=min_degree, loops=loops):
            h = h.canonical_label()
            g = DirectedGraph(num_vertices, list(h.edges(labels=False)))
            g.canonicalize_edges()
            if has_odd_automorphism is None or directed_graph_has_odd_automorphism(g) == has_odd_automorphism:
                yield g

def directed_graph_to_encoding(g):
    n = len(g)
    edges = g.edges()
    G = DiGraph([list(range(n)), edges])
    return '&' + G.dig6_string()

def directed_graph_from_encoding(digraph6_string):
    G = DiGraph(digraph6_string[1:])
    return DirectedGraph(len(G.vertices()), list(G.edges(labels=False)))
