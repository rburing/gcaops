import itertools
import subprocess
import os
import sage.all # make SageMath work when called from Python
from sage.env import SAGE_NAUTY_BINS_PREFIX
from sage.graphs.digraph import DiGraph
from gcaops.graph.formality_graph import FormalityGraph
from .permutation import selection_sort

def nauty_generate_formality(num_ground_vertices, num_aerial_vertices, num_undirected_edges, num_directed_edges,
        connected=None, max_out_degree=None, num_verts_of_max_out_degree=None, loops=None):
    num_vertices = num_ground_vertices + num_aerial_vertices
    geng_args = [str(num_vertices), "{}:{}".format(num_undirected_edges, num_undirected_edges)]
    if connected is not None and connected:
        geng_args.append("-c")
    directg_args = ['-e{}:{}'.format(num_directed_edges, num_directed_edges)]
    if loops is not None and not loops:
        directg_args.append('-o')
    pickg_args = []
    if num_ground_vertices > 0: # NOTE: avoid edge case of pickg seemingly not working
        pickg_args = ['-d0', '-m{}:'.format(num_ground_vertices)] # NOTE: can have more vertices of out-degree 0, namely aerial vertices
    if max_out_degree:
        pickg_args.append('-D{}'.format(max_out_degree))
    if num_verts_of_max_out_degree:
        pickg_args.append('-M{}'.format(num_verts_of_max_out_degree))
    FNULL = open(os.devnull, 'w')
    geng = subprocess.Popen((SAGE_NAUTY_BINS_PREFIX + 'geng', *geng_args), stdout=subprocess.PIPE, stderr=FNULL)
    directg = subprocess.Popen((SAGE_NAUTY_BINS_PREFIX + 'directg', *directg_args), stdin=geng.stdout, stdout=subprocess.PIPE, stderr=FNULL)
    pickg = subprocess.Popen((SAGE_NAUTY_BINS_PREFIX + 'pickg', *pickg_args), stdin=directg.stdout, stdout=subprocess.PIPE, stderr=FNULL)
    for line in pickg.stdout:
        digraph6_string = line.decode('ascii').rstrip()
        yield DiGraph(digraph6_string[1:])

def formality_graph_generate(num_ground_vertices, num_aerial_vertices, num_edges,
        connected=None, max_out_degree=None, num_verts_of_max_out_degree=None, sorted_out_degrees=None,
        max_aerial_in_degree=None,
        loops=None, prime=None,
        has_odd_automorphism=None, positive_differential_order=None, mod_ground_permutations=False):
    if loops is None or loops:
        max_loop_order = num_edges // 2 # NOTE: can have at most this many loops, while still attaining num_edges
    else:
        max_loop_order = 0
    num_vertices = num_ground_vertices + num_aerial_vertices
    partition = [[v] for v in range(num_ground_vertices)] + [list(range(num_ground_vertices, num_vertices))]
    for loop_order in range(max_loop_order + 1):
        for h in nauty_generate_formality(num_ground_vertices, num_aerial_vertices, num_edges - loop_order, num_edges,
                                          connected=connected, max_out_degree=max_out_degree,
                                          num_verts_of_max_out_degree=num_verts_of_max_out_degree, loops=loops):
            if sorted_out_degrees is not None and tuple(sorted(h.out_degree_sequence())) != sorted_out_degrees:
                continue
            # Choose sinks
            possible_sinks = [v for v in h if h.out_degree(v) == 0]
            # TODO: instead of all combinations, mod out by automorphisms
            seen = []
            for sinks in itertools.combinations(possible_sinks, num_ground_vertices):
                non_sinks = tuple(v for v in h if not v in sinks)
                if max_aerial_in_degree is not None and max(h.in_degree(v) for v in non_sinks) > max_aerial_in_degree:
                    continue
                # Relabel sinks to 0, 1, ...
                relabeling = dict(zip(sinks + non_sinks, range(num_vertices)))
                k = h.relabel(relabeling, inplace=False)
                if mod_ground_permutations:
                    ground_permutations = [list(range(num_ground_vertices))]
                else:
                    ground_permutations = itertools.permutations(range(num_ground_vertices))
                for sigma in ground_permutations:
                    hh = k.relabel(dict(zip(range(num_ground_vertices), sigma)), inplace=False)
                    hh = hh.canonical_label(partition=partition)
                    g = FormalityGraph(num_ground_vertices, num_aerial_vertices, list(hh.edges(labels=False, sort=True)))
                    g.canonicalize_edges()
                    if g in seen:
                        continue
                    else:
                        seen.append(g)
                    if positive_differential_order is not None and positive_differential_order != (not 0 in g.differential_orders()):
                        continue
                    if has_odd_automorphism is not None and formality_graph_has_odd_automorphism(g) != has_odd_automorphism:
                        continue
                    if prime is not None and formality_graph_is_prime(g) != prime:
                        continue
                    yield g

def formality_graph_canonicalize(g):
    n = len(g)
    edges = g.edges()
    partition = [[v] for v in range(g.num_ground_vertices())] + [list(range(g.num_ground_vertices(), g.num_ground_vertices() + g.num_aerial_vertices()))]
    H = DiGraph([list(range(n)), edges], loops=True)
    if len(H.edges()) != len(g.edges()):
        raise ValueError("don't know how to canonicalize graph with double edges")
    G, sigma = H.canonical_label(partition=partition, certificate=True)
    new_edges = list(G.edges(labels=False, sort=True))
    edge_permutation = [tuple([sigma[edge[0]],sigma[edge[1]]]) for edge in edges]
    index_permutation = [new_edges.index(e) for e in edge_permutation]
    undo_canonicalize = [0]*n
    for k, v in sigma.items():
        undo_canonicalize[v] = k
    return FormalityGraph(g.num_ground_vertices(), g.num_aerial_vertices(), list(new_edges)), undo_canonicalize, selection_sort(index_permutation)

def formality_graph_has_odd_automorphism(g):
    n = len(g)
    edges = g.edges()
    partition = [[v] for v in range(g.num_ground_vertices())] + [list(range(g.num_ground_vertices(), g.num_ground_vertices() + g.num_aerial_vertices()))]
    G = DiGraph([list(range(n)), edges])
    for sigma in G.automorphism_group(partition=partition).gens(): # NOTE: it suffices to check generators
        edge_permutation = [tuple([sigma(edge[0]),sigma(edge[1])]) for edge in edges]
        index_permutation = [edges.index(e) for e in edge_permutation]
        if selection_sort(index_permutation) == -1:
            return True
    return False

def formality_graph_is_prime(g):
    n = len(g)
    G = DiGraph([list(range(n)), g.edges()])
    G.delete_vertices(list(range(g.num_ground_vertices())))
    return G.is_connected()

def formality_graph_to_encoding(g):
    n = len(g)
    edges = g.edges()
    G = DiGraph([list(range(n)), edges])
    return (g.num_ground_vertices(), g.num_aerial_vertices(), G.dig6_string()) # TODO: do better?

def formality_graph_from_encoding(encoding):
    num_ground, num_aerial, dig6_string = encoding
    G = DiGraph(dig6_string)
    return FormalityGraph(num_ground, num_aerial, list(G.edges(labels=False, sort=True)))
