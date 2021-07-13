from .graph_file import UndirectedGraphFileView, DirectedGraphFileView, UndirectedToDirectedGraphFileView, FormalityGraphFileView
from abc import ABC, abstractmethod
from util.undirected_graph_sage import undirected_graph_canonicalize, undirected_graph_generate
from util.directed_graph_sage import directed_graph_canonicalize, directed_graph_generate_from_undirected, undirected_to_directed_graph_coefficient
from util.formality_graph_sage import formality_graph_canonicalize, formality_graph_generate
import os

GRAPH_CACHE_DIR = None # e.g. '/home/rburing/src/gcaops/.graph_cache'

# TODO: Reuse graphs to save memory.

def options_to_filename(vertices, edges, directed=False, connected=False, biconnected=False, min_degree=0, loops=None, has_odd_automorphism=None):
    filename = '{}_{}_{}'.format('d' if directed else 'u', vertices, edges)
    if connected:
        filename += '_c'
    if biconnected:
        filename += '_bc'
    if min_degree != 0:
        filename += '_m{}'.format(min_degree)
    if loops is not None:
        filename += '_loops' if loops else '_noloops'
    if has_odd_automorphism is not None:
        filename += '_odd' if has_odd_automorphism else '_even'
    filename += '.db'
    return filename

def formality_options_to_filename(num_ground_vertices, num_aerial_vertices, num_edges, connected=None,
                                  max_out_degree=None, num_verts_of_max_out_degree=None, loops=None, prime=None,
                                  has_odd_automorphism=None, positive_differential_order=None, mod_ground_permutations=False):
    filename = 'f_{}_{}_{}'.format(num_ground_vertices, num_aerial_vertices, num_edges)
    if connected is not None:
        filename += '_c' if connected else '_nc'
    if max_out_degree is not None:
        filename += '_D{}'.format(max_out_degree)
    if num_verts_of_max_out_degree is not None:
        filename += '_M{}'.format(num_verts_of_max_out_degree)
    if loops is not None:
        filename += '_loops' if loops else '_noloops'
    if prime is not None:
        filename += '_prime' if prime else '_composite'
    if has_odd_automorphism is not None:
        filename += '_odd' if has_odd_automorphism else '_even'
    if positive_differential_order is not None:
        filename += '_pdo' if positive_differential_order else '_npdo'
    if mod_ground_permutations:
        filename += '_modground'
    filename += '.db'
    return filename

class GraphCache(ABC):
    file_view = None.__class__
    cache_keys = []
    cache = {}

    @abstractmethod
    def canonicalize_graph(self, graph):
        pass

    @abstractmethod
    def graphs(self, grading, **options):
        pass

class UndirectedGraphCache(GraphCache):
    cache_keys = ['connected', 'biconnected', 'min_degree', 'has_odd_automorphism']
    file_view = UndirectedGraphFileView

    def canonicalize_graph(self, graph):
        return undirected_graph_canonicalize(graph)

    def _add_graphs(self, result, bi_grading, **options):
        if len(result) != 0:
            return # assume we've already done this (valid assumption because commit is done *only after* adding *all* the graphs)

        num_vertices, num_edges = bi_grading
        del options['directed']
        for g in undirected_graph_generate(num_vertices, num_edges, **options):
            result.append(g)

    def graphs(self, bi_grading, connected=False, biconnected=False, min_degree=0, has_odd_automorphism=True):
        options = {'directed': False, 'connected': connected, 'biconnected': biconnected, 'min_degree': min_degree, 'has_odd_automorphism': has_odd_automorphism}
        num_vertices, num_edges = bi_grading
        cache_key = bi_grading + tuple(options[k] for k in self.cache_keys)
        if cache_key in self.cache:
            return self.cache[cache_key]
        if GRAPH_CACHE_DIR is not None:
            if not os.path.isdir(GRAPH_CACHE_DIR):
                os.makedirs(GRAPH_CACHE_DIR)
            filename = os.path.join(GRAPH_CACHE_DIR, options_to_filename(num_vertices, num_edges, **options))
            my_graphs = self.file_view(filename, num_vertices, num_edges)
        else:
            my_graphs = list()
        self._add_graphs(my_graphs, bi_grading, **options)
        if GRAPH_CACHE_DIR is not None:
            my_graphs.commit()
        self.cache[cache_key] = my_graphs
        return my_graphs

class DirectedGraphCache(GraphCache):
    cache_keys = ['connected', 'biconnected', 'min_degree', 'loops', 'has_odd_automorphism']
    file_view = DirectedGraphFileView

    def __init__(self, undirected_graph_cache):
        self._undirected_graph_cache = undirected_graph_cache
        self._undirected_to_directed = {}

    def canonicalize_graph(self, graph):
        return directed_graph_canonicalize(graph)

    def _add_graphs(self, result, orientation_data, bi_grading, **options):
        if len(result) != 0:
            return # assume we've already done this (valid assumption because commit is done *only after* adding *all* the graphs)

        num_vertices, num_edges = bi_grading

        # NOTE: for an even directed graph (having no odd automorphism), the underlying undirected graph may still be odd (have an odd automorphism)
        underlying_undirected_options = {'connected' : options['connected'], 'biconnected' : options['biconnected'], 'min_degree' : options['min_degree'], 'has_odd_automorphism' : None}

        undirected_options = {}
        undirected_options.update(underlying_undirected_options)
        undirected_options['has_odd_automorphism'] = options['has_odd_automorphism']

        if options['loops'] is None or options['loops']:
            max_loop_order = num_edges // 2 # NOTE: can have at most this many loops, while still attaining num_edges
        else:
            max_loop_order = 0

        for loop_order in range(max_loop_order + 1):
            loopless = loop_order == 0
            h_idx = 0
            # TODO: use iterator over encodings?
            for g in self._undirected_graph_cache.graphs((num_vertices, num_edges - loop_order), **underlying_undirected_options):
                # TODO: arrange to get the index of the undirected graph (for orientation data) more efficiently
                try:
                    if loopless: # NOTE: only graphs without loops are in the image of the orientation map
                        g_idx = self._undirected_graph_cache.graphs((num_vertices, num_edges), **undirected_options).index(g)
                    else:
                        g_idx = None
                except ValueError:
                    g_idx = None
                for h in directed_graph_generate_from_undirected(g, num_edges, loops=options['loops'], has_odd_automorphism=options['has_odd_automorphism']):
                    result.append(h)
                    if g_idx is not None:
                        c = undirected_to_directed_graph_coefficient(g, h)
                        orientation_data.append((g_idx, h_idx, c))
                    h_idx += 1

    def graphs(self, bi_grading, connected=False, biconnected=False, min_degree=0, loops=True, has_odd_automorphism=True):
        options = {'directed': True, 'connected': connected, 'biconnected': biconnected, 'min_degree': min_degree, 'loops': loops, 'has_odd_automorphism': has_odd_automorphism}
        num_vertices, num_edges = bi_grading
        cache_key = bi_grading + tuple(options[k] for k in self.cache_keys)
        if cache_key in self.cache:
            return self.cache[cache_key]
        if GRAPH_CACHE_DIR is not None:
            if not os.path.isdir(GRAPH_CACHE_DIR):
                os.makedirs(GRAPH_CACHE_DIR)
            basename = options_to_filename(num_vertices, num_edges, **options)
            filename = os.path.join(GRAPH_CACHE_DIR, basename)
            orientation_filename = os.path.join(GRAPH_CACHE_DIR, 'u_to_' + basename)
            orientation_data = UndirectedToDirectedGraphFileView(orientation_filename)
            my_graphs = self.file_view(filename, num_vertices, num_edges)
        else:
            orientation_data = list()
            my_graphs = list()
        self._add_graphs(my_graphs, orientation_data, bi_grading, **options)
        if GRAPH_CACHE_DIR is not None:
            my_graphs.commit()
            orientation_data.commit()
        self.cache[cache_key] = my_graphs
        self._undirected_to_directed[cache_key] = orientation_data
        return my_graphs

    def _undirected_to_directed_coeffs(self, bi_grading, undirected_graph_idx, **options):
        cache_key = bi_grading + tuple(options[k] for k in self.cache_keys)
        num_vertices, num_edges = bi_grading
        self.graphs(bi_grading, **options) # NOTE: this makes sure the required orientation data has been generated/loaded
        orientation_data = self._undirected_to_directed[cache_key]
        if GRAPH_CACHE_DIR is not None:
            yield from orientation_data.undirected_to_directed_coeffs(undirected_graph_idx)
        else:
            for row in orientation_data:
                if row[0] == undirected_graph_idx:
                    yield (row[1], row[2])

class FormalityGraphCache(GraphCache):
    cache_keys = ['connected', 'max_out_degree', 'num_verts_of_max_out_degree', 'loops', 'prime', 'has_odd_automorphism', 'positive_differential_order', 'mod_ground_permutations']
    file_view = FormalityGraphFileView

    def canonicalize_graph(self, graph):
        return formality_graph_canonicalize(graph)

    def _add_graphs(self, result, tri_grading, **options):
        if len(result) != 0:
            return # assume we've already done this (valid assumption because commit is done *only after* adding *all* the graphs)

        num_ground_vertices, num_aerial_vertices, num_edges = tri_grading
        for g in formality_graph_generate(num_ground_vertices, num_aerial_vertices, num_edges, **options):
            result.append(g)

    def graphs(self, tri_grading, connected=None, max_out_degree=None, num_verts_of_max_out_degree=None, loops=None, prime=None,
               has_odd_automorphism=None, positive_differential_order=None, mod_ground_permutations=False):
        options = {'connected': connected, 'max_out_degree' : max_out_degree, 'num_verts_of_max_out_degree' : num_verts_of_max_out_degree,
                   'loops': loops, 'prime': prime, 'has_odd_automorphism': has_odd_automorphism, 'positive_differential_order': positive_differential_order,
                   'mod_ground_permutations' : mod_ground_permutations}
        num_ground_vertices, num_aerial_vertices, num_edges = tri_grading
        cache_key = tri_grading + tuple(options[k] for k in self.cache_keys)
        if cache_key in self.cache:
            return self.cache[cache_key]
        if GRAPH_CACHE_DIR is not None:
            if not os.path.isdir(GRAPH_CACHE_DIR):
                os.makedirs(GRAPH_CACHE_DIR)
            filename = os.path.join(GRAPH_CACHE_DIR, formality_options_to_filename(num_ground_vertices, num_aerial_vertices, num_edges, **options))
            my_graphs = self.file_view(filename, num_ground_vertices, num_aerial_vertices, num_edges)
        else:
            my_graphs = list()
        self._add_graphs(my_graphs, tri_grading, **options)
        if GRAPH_CACHE_DIR is not None:
            my_graphs.commit()
        self.cache[cache_key] = my_graphs
        return my_graphs

undirected_graph_cache = UndirectedGraphCache()
directed_graph_cache = DirectedGraphCache(undirected_graph_cache)
formality_graph_cache = FormalityGraphCache()
