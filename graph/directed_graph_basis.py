from graph.directed_graph import DirectedGraph
from graph.graph_basis import GraphBasis
from graph.graph_cache import directed_graph_cache
from util.misc import keydefaultdict
from functools import partial

class DirectedGraphBasis(GraphBasis):
    """
    Basis of a module spanned by directed graphs.
    
    A basis consists of keys ``(v,e,index,...)`` where ``(v,e,index)`` identifies the isomorphism class of the graph.
    """
    graph_class = DirectedGraph
    grading_size = 2

class DirectedGraphComplexBasis(DirectedGraphBasis):
    """
    Basis consisting of representatives of isomorphism classes of directed graphs with no automorphisms that induce an odd permutation on edges
    """
    def __init__(self, connected=None, biconnected=None, min_degree=0, loops=True):
        """
        Initialize this basis.
        """
        if not min_degree in [0, 3]:
            raise ValueError('min_degree can only be 0 or 3')
        self._connected = connected
        self._biconnected = biconnected
        self._min_degree = min_degree
        self._loops = loops
        self._graphs = keydefaultdict(partial(directed_graph_cache.graphs, connected=connected, biconnected=biconnected, min_degree=min_degree, loops=loops, has_odd_automorphism=False))

    def graph_to_key(self, graph):
        """
        Return a tuple consisting of the key in this basis and the sign factor such that ``graph`` equals the sign times the graph identified by the key.

        INPUT:

        - ``graph`` -- a DirectedGraph

        OUTPUT:

        Either ``(None, 1)`` if the input ``graph`` is not in the span of the basis, or a tuple consisting of a key and a sign, where a key is a tuple consisting of the number of vertices, the number of edges, and the index of the graph in the list.
        """
        g, _, sign = directed_graph_cache.canonicalize_graph(graph)
        v, e = len(g), len(g.edges())
        try:
            index = self._graphs[v,e].index(g)
            return (v,e,index), sign
        except ValueError:
            return None, 1

    def key_to_graph(self, key):
        """
        Return a tuple consisting of a DirectedGraph and the sign factor such that the sign times the graph equals the graph identified by the key.

        INPUT:

        - ``key`` -- a key in this basis

        OUTPUT:

        Either ``(None, 1)`` if the input ``key`` is not in the basis, or a tuple consisting of a DirectedGraph and a sign which is always +1.
        """
        v, e, index = key 
        try:
            return self._graphs[v,e][index], 1
        except IndexError:
            return None, 1

    def __repr__(self):
        """
        Return a string representation of this basis.
        """
        filters = []
        if self._connected:
            filters.append('connected')
        if self._biconnected:
            filters.append('biconnected')
        if self._min_degree != 0:
            filters.append('of degree at least {}'.format(self._min_degree))
        if not self._loops:
            filters.append('without loops')
        if filters:
            filters_str = ' ({})'.format(', '.join(filters))
        else:
            filters_str = ''
        return 'Basis consisting of representatives of isomorphism classes of directed graphs{} with no automorphisms that induce an odd permutation on edges'.format(filters_str)

    def graph_properties(self):
        """
        Return a dictionary containing the properties of the graphs in this basis.
        """
        return {'connected' : self._connected, 'biconnected' : self._biconnected, 'min_degree' : self._min_degree, 'loops' : self._loops, 'has_odd_automorphism' : False}

    def graphs(self, vertices, edges):
        """
        Return the list of graphs in this basis with the given amount of ``vertices`` and ``edges``.
        """
        return self._graphs[vertices, edges]

    def cardinality(self, vertices, edges):
        """
        Return the number of graphs in this basis with the given amount of ``vertices`` and ``edges``.
        """
        return len(self._graphs[vertices,edges])

    def _undirected_to_directed_coeffs(self, bi_grading, undirected_graph_idx):
        """
        Return an iterator over tuples of the form ``(directed_graph_idx, coefficient)``, where the index is an index in this basis.

        These tuples correspond to terms in the orientation of the respective undirected graph.

        ASSUMPTIONS:

        Assumes that ``undirected_graph_idx`` refers to an index in a basis with the same options as ``self`` (except for the ``loops`` option).
        """
        options = self.graph_properties()
        yield from directed_graph_cache._undirected_to_directed_coeffs(bi_grading, undirected_graph_idx, **options)

class DirectedGraphOperadBasis(DirectedGraphBasis):
    """
    Basis consisting of labeled directed graphs with no automorphisms that induce an odd permutation on edges
    """
    def __init__(self):
        """
        Initialize this basis.
        """
        self._graphs = keydefaultdict(partial(directed_graph_cache.graphs, has_odd_automorphism=False))

    def graph_to_key(self, graph):
        """
        Return a tuple consisting of the key in this basis and the sign factor such that ``graph`` equals the sign times the graph identified by the key.

        INPUT:

        - ``graph`` -- a DirectedGraph

        OUTPUT:

        Either ``(None, 1)`` if the input ``graph`` is not in the span of the basis, or a tuple consisting of a key and a sign, where a key is a tuple consisting of the number of vertices, the number of edges, the index of the graph in the list, followed by a permutation of vertices.
        """
        g, undo_canonicalize, sign = directed_graph_cache.canonicalize_graph(graph)
        v, e = len(g), len(g.edges())
        try:
            index = self._graphs[v,e].index(g)
            return (v,e,index) + tuple(undo_canonicalize), sign
        except ValueError:
            return None, 1

    def key_to_graph(self, key):
        """
        Return a tuple consisting of a DirectedGraph and the sign factor such that the sign times the graph equals the graph identified by the key.

        INPUT:

        - ``key`` -- a key in this basis

        OUTPUT:

        Either ``(None, 1)`` if the input ``key`` is not in the basis, or a tuple consisting of a DirectedGraph and a sign.
        """
        v, e, index = key[:3]
        undo_canonicalize = key[3:]
        try:
            G = self._graphs[v,e][index]
            g = G.relabeled(undo_canonicalize)
            sign = g.canonicalize_edges()
            return g, sign
        except IndexError:
            return None, 1

    def __repr__(self):
        """
        Return a string representation of this basis.
        """
        return 'Basis consisting of labeled directed graphs with no automorphisms that induce an odd permutation on edges'

    def graph_properties(self):
        """
        Return a dictionary containing the properties of the graphs in this basis.
        """
        return {'has_odd_automorphism' : False}
