from .graph_file import UndirectedGraphFileView, DirectedGraphFileView
from abc import ABC, abstractmethod
from util.undirected_graph_sage import undirected_graph_canonicalize, undirected_graph_generate
from util.directed_graph_sage import directed_graph_canonicalize, directed_graph_generate_from_undirected, undirected_to_directed_graph_coefficient
import os
import sqlite3

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

class GraphCache(ABC):
    cache_keys = []
    cache = {}

    @abstractmethod
    def canonicalize_graph(self, graph):
        pass

    @abstractmethod
    def _add_graphs(self, result, bi_grading, **options):
        pass

    def graphs(self, bi_grading, **options):
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

class UndirectedGraphCache(GraphCache):
    cache_keys = ['connected', 'biconnected', 'min_degree', 'has_odd_automorphism']
    file_view = UndirectedGraphFileView

    def canonicalize_graph(self, graph):
        return undirected_graph_canonicalize(graph)

    def _add_graphs(self, result, bi_grading, **options):
        num_vertices, num_edges = bi_grading
        del options['directed']
        if len(result) == 0:
            for g in undirected_graph_generate(num_vertices, num_edges, **options):
                result.append(g)

    def graphs(self, bi_grading, connected=False, biconnected=False, min_degree=0, has_odd_automorphism=True):
        options = {'directed': False, 'connected': connected, 'biconnected': biconnected, 'min_degree': min_degree, 'has_odd_automorphism': has_odd_automorphism}
        return super().graphs(bi_grading, **options)

class DirectedGraphCache(GraphCache):
    cache_keys = ['connected', 'biconnected', 'min_degree', 'loops', 'has_odd_automorphism']
    file_view = DirectedGraphFileView

    def __init__(self, undirected_graph_cache):
        self._undirected_graph_cache = undirected_graph_cache

    def canonicalize_graph(self, graph):
        return directed_graph_canonicalize(graph)

    def _add_graphs(self, result, bi_grading, **options):
        num_vertices, num_edges = bi_grading
        undirected_options = {'connected' : options['connected'], 'biconnected' : options['biconnected'], 'min_degree' : options['min_degree'], 'has_odd_automorphism' : options['has_odd_automorphism']}
        if len(result) == 0:
            if options['loops'] is None or options['loops']:
                max_loop_order = num_edges // 2 # NOTE: can have at most this many loops, while still attaining num_edges
            else:
                max_loop_order = 0

            if GRAPH_CACHE_DIR is not None:
                basename = options_to_filename(num_vertices, num_edges, **options)
                orientation_filename = os.path.join(GRAPH_CACHE_DIR, 'u_to_' + basename)
                orientation_con = sqlite3.connect(orientation_filename)
                orientation_cur = orientation_con.cursor()
                orientation_cur.execute('DROP TABLE IF EXISTS undirected_to_directed')
                orientation_cur.execute('CREATE TABLE undirected_to_directed (undirected_graph_id INTEGER, directed_graph_id INTEGER, coefficient INTEGER)')
                orientation_cur.execute('CREATE INDEX index_undirected ON undirected_to_directed(undirected_graph_id)')
                orientation_con.commit()

            for loop_order in range(max_loop_order + 1):
                # TODO: use iterator over encodings?
                for (g_idx, g) in enumerate(self._undirected_graph_cache.graphs((num_vertices, num_edges - loop_order), **undirected_options)):
                    for (h_idx,h) in enumerate(directed_graph_generate_from_undirected(g, num_edges, loops=options['loops'], has_odd_automorphism=options['has_odd_automorphism'])):
                        result.append(h)
                        if GRAPH_CACHE_DIR is not None and loop_order == 0: # NOTE: only graphs without loops are in the image of the orientation map
                            c = undirected_to_directed_graph_coefficient(g, h)
                            orientation_cur.execute('INSERT INTO undirected_to_directed (undirected_graph_id, directed_graph_id, coefficient) VALUES (?, ?, ?)', (1 + g_idx, 1 + h_idx, c))
            if GRAPH_CACHE_DIR is not None:
                orientation_con.commit()

    def graphs(self, bi_grading, connected=False, biconnected=False, min_degree=0, loops=True, has_odd_automorphism=True):
        options = {'directed': True, 'connected': connected, 'biconnected': biconnected, 'min_degree': min_degree, 'loops': loops, 'has_odd_automorphism': has_odd_automorphism}
        return super().graphs(bi_grading, **options)

undirected_graph_cache = UndirectedGraphCache()
directed_graph_cache = DirectedGraphCache(undirected_graph_cache)
