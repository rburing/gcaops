from graph.directed_graph import DirectedGraph
from util.undirected_graph_sage import undirected_graph_to_encoding
from util.permutation import selection_sort
from itertools import combinations
import sage.all # make SageMath work when called from Python
from sage.graphs.digraph import DiGraph
from sage.graphs.graph import Graph
import subprocess
import os

NAUTY_PREFIX = '' # e.g. '/home/rburing/src/nauty27r1/'

# TODO: use a single directg process for multiple inputs

def nauty_generate_directed_from_undirected(undirected_graph, num_directed_edges, loops=True):
    num_vertices = len(undirected_graph)
    directg_args = ['-T', '-e{}:{}'.format(num_directed_edges, num_directed_edges)]
    if not loops:
        directg_args.append('-o')
    FNULL = open(os.devnull, 'w')
    directg = subprocess.Popen((NAUTY_PREFIX + 'directg', *directg_args), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=FNULL, encoding='ascii')
    outs, errs = directg.communicate(undirected_graph_to_encoding(undirected_graph) + '\n')
    for line in outs.split('\n'):
        graph_encoding = line.rstrip()
        if graph_encoding == '':
            continue
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

def directed_graph_generate_from_undirected(undirected_graph, num_directed_edges, loops=True, has_odd_automorphism=None):
    num_vertices = len(undirected_graph)
    for h in nauty_generate_directed_from_undirected(undirected_graph, num_directed_edges, loops=loops):
        h = h.canonical_label()
        g = DirectedGraph(num_vertices, list(h.edges(labels=False)))
        g.canonicalize_edges()
        if has_odd_automorphism is None or directed_graph_has_odd_automorphism(g) == has_odd_automorphism:
            yield g

def undirected_to_directed_graph_coefficient(undirected_graph, directed_graph):
    g = Graph(undirected_graph.edges())
    h = Graph(directed_graph.edges())
    are_isomorphic, sigma = g.is_isomorphic(h, certificate=True)
    assert are_isomorphic
    edges = directed_graph.edges()
    edge_permutation = [edges.index((sigma[a], sigma[b])) if (sigma[a], sigma[b]) in edges else edges.index((sigma[b], sigma[a])) for (a,b) in undirected_graph.edges()]
    sign = selection_sort(edge_permutation)
    multiplicity = len(g.automorphism_group()) // len(DiGraph(directed_graph.edges()).automorphism_group())
    return sign * multiplicity

def directed_graph_to_encoding(g):
    n = len(g)
    edges = g.edges()
    G = DiGraph([list(range(n)), edges])
    return '&' + G.dig6_string()

def directed_graph_from_encoding(digraph6_string):
    G = DiGraph(digraph6_string[1:])
    return DirectedGraph(len(G.vertices()), list(G.edges(labels=False)))
