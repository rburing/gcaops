from abc import ABC, abstractmethod
from util.undirected_graph_sage import undirected_graph_canonicalize, undirected_graph_generate
from util.directed_graph_sage import directed_graph_canonicalize, directed_graph_generate
import os
import pickle

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
    if loops:
        filename += '_loops'
    if has_odd_automorphism is not None:
        filename += '_odd' if has_odd_automorphism else '_even'
    filename += '.pickle'
    return filename

class GraphCache(ABC):
    cache_keys = []
    cache = {}

    @abstractmethod
    def canonicalize_graph(self, graph):
        pass

    @abstractmethod
    def graphs(self, bi_grading, **options):
        pass

    @abstractmethod
    def _generate_graphs(self, bi_grading, **options):
        pass

    def graphs(self, bi_grading, **options):
        vertices, edges = bi_grading
        cache_key = bi_grading + tuple(options[k] for k in self.cache_keys)
        if cache_key in self.cache:
            return self.cache[cache_key]
        if GRAPH_CACHE_DIR is not None:
            if not os.path.isdir(GRAPH_CACHE_DIR):
                os.makedirs(GRAPH_CACHE_DIR)
            filename = os.path.join(GRAPH_CACHE_DIR, options_to_filename(vertices, edges, **options))
            if os.path.isfile(filename):
                with open(filename, 'rb') as f:
                    my_graphs = pickle.load(f)
            else:
                my_graphs = list(self._generate_graphs(bi_grading, **options))
                with open(filename, 'wb') as f:
                    pickle.dump(my_graphs, f, pickle.HIGHEST_PROTOCOL)
        else:
            my_graphs = list(self._generate_graphs(bi_grading, **options))
        self.cache[cache_key] = my_graphs
        return my_graphs

class UndirectedGraphCache(GraphCache):
    cache_keys = ['directed', 'connected', 'biconnected', 'min_degree', 'has_odd_automorphism']

    def canonicalize_graph(self, graph):
        return undirected_graph_canonicalize(graph)

    def _generate_graphs(self, bi_grading, **options):
        vertices, edges = bi_grading
        del options['directed']
        for g in undirected_graph_generate(vertices, edges, **options):
            yield g

    def graphs(self, bi_grading, connected=False, biconnected=False, min_degree=0, has_odd_automorphism=True):
        options = {'directed': False, 'connected': connected, 'biconnected': biconnected, 'min_degree': min_degree, 'has_odd_automorphism': has_odd_automorphism}
        return super().graphs(bi_grading, **options)

class DirectedGraphCache(GraphCache):
    cache_keys = ['directed', 'connected', 'biconnected', 'min_degree', 'loops', 'has_odd_automorphism']

    def canonicalize_graph(self, graph):
        return directed_graph_canonicalize(graph)

    def _generate_graphs(self, bi_grading, **options):
        vertices, edges = bi_grading
        del options['directed']
        for g in directed_graph_generate(vertices, edges, **options):
            yield g

    def graphs(self, bi_grading, connected=False, biconnected=False, min_degree=0, loops=True, has_odd_automorphism=True):
        options = {'directed': True, 'connected': connected, 'biconnected': biconnected, 'min_degree': min_degree, 'loops': loops, 'has_odd_automorphism': has_odd_automorphism}
        return super().graphs(bi_grading, **options)

undirected_graph_cache = UndirectedGraphCache()
directed_graph_cache = DirectedGraphCache()
