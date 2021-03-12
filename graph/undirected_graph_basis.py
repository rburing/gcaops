from util.misc import keydefaultdict
from util.permutation import selection_sort
from util.undirected_graph_sage import undirected_graph_canonicalize, undirected_graph_generate, undirected_graph_has_odd_automorphism
from abc import ABC, abstractmethod

class UndirectedGraphBasis(ABC):
    """
    Basis of a module of undirected graphs.
    
    A basis consists of keys ``(v,e,index,...)`` where ``(v,e,index)`` identifies the isomorphism class of the graph.
    """
    @abstractmethod
    def graph_to_key(self, graph):
        """
        Return a tuple consisting of the key in ``self`` and the sign factor such that ``graph`` equals the sign times the graph identified by the key.

        INPUT:

        - ``graph`` -- an UndirectedGraph
        """
        pass

    @abstractmethod
    def key_to_graph(self, key):
        """
        Return a tuple consisting of an UndirectedGraph and the sign factor such that the sign times the graph equals the graph identified by the key.

        INPUT:

        - ``key`` -- a key in ``self``
        """
        pass

    def __repr__(self):
        """
        Return a string representation of ``self``.
        """
        return 'Basis consisting of undirected graphs'

class UndirectedGraphComplexBasis(UndirectedGraphBasis):
    """
    Basis consisting of representatives of isomorphism classes of undirected graphs with no automorphisms that induce an odd permutation on edges
    """
    def __init__(self, connected=None, biconnected=None, min_degree=0):
        """
        Initialize ``self``.
        """
        self._connected = connected
        self._biconnected = biconnected
        self._min_degree = min_degree
        self._graphs = keydefaultdict(lambda key: list(filter(lambda g: not undirected_graph_has_odd_automorphism(g), undirected_graph_generate(*key, connected=connected, biconnected=biconnected, min_degree=min_degree))))

    def graph_to_key(self, graph):
        """
        Return a tuple consisting of the key in ``self`` and the sign factor such that ``graph`` equals the sign times the graph identified by the key.

        INPUT:

        - ``graph`` -- an UndirectedGraph

        OUTPUT:

        Either ``(None, 1)`` if the input ``graph`` is not in the span of the basis, or a tuple consisting of a key and a sign, where a key is a tuple consisting of the number of vertices, the number of edges, and the index of the graph in the list.
        """
        g, _, sign = undirected_graph_canonicalize(graph)
        v, e = len(g), len(g.edges())
        try:
            index = self._graphs[v,e].index(g)
            return (v,e,index), sign
        except ValueError:
            return None, 1

    def key_to_graph(self, key):
        """
        Return a tuple consisting of an UndirectedGraph and the sign factor such that the sign times the graph equals the graph identified by the key.

        INPUT:

        - ``key`` -- a key in ``self``

        OUTPUT:

        Either ``(None, 1)`` if the input ``key`` is not in the basis, or a tuple consisting of an UndirectedGraph and a sign which is always +1.
        """
        v, e, index = key 
        try:
            return self._graphs[v,e][index], 1
        except IndexError:
            return None, 1

    def __repr__(self):
        """
        Return a string representation of ``self``.
        """
        filters = []
        if self._connected:
            filters.append('connected')
        if self._biconnected:
            filters.append('biconnected')
        if self._min_degree != 0:
            filters.append('of degree at least {}'.format(self._min_degree))
        if filters:
            filters_str = ' ({})'.format(', '.join(filters))
        else:
            filters_str = ''
        return 'Basis consisting of representatives of isomorphism classes of undirected graphs{} with no automorphisms that induce an odd permutation on edges'.format(filters_str)

class UndirectedGraphOperadBasis(UndirectedGraphBasis):
    """
    Basis consisting of representatives of isomorphism classes of labeled undirected graphs with no automorphisms that induce an odd permutation on edges
    """
    def __init__(self):
        """
        Initialize ``self``.
        """
        self._graphs = keydefaultdict(lambda key: list(filter(lambda g: not undirected_graph_has_odd_automorphism(g), undirected_graph_generate(*key))))

    def graph_to_key(self, graph):
        """
        Return a tuple consisting of the key in ``self`` and the sign factor such that ``graph`` equals the sign times the graph identified by the key.

        INPUT:

        - ``graph`` -- an UndirectedGraph

        OUTPUT:

        Either ``(None, 1)`` if the input ``graph`` is not in the span of the basis, or a tuple consisting of a key and a sign, where a key is a tuple consisting of the number of vertices, the number of edges, the index of the graph in the list, followed by a permutation of vertices.
        """
        g, undo_canonicalize, sign = undirected_graph_canonicalize(graph)
        v, e = len(g), len(g.edges())
        try:
            index = self._graphs[v,e].index(g)
            return (v,e,index) + tuple(undo_canonicalize), sign
        except ValueError:
            return None, 1

    def key_to_graph(self, key):
        """
        Return a tuple consisting of an UndirectedGraph and the sign factor such that the sign times the graph equals the graph identified by the key.

        INPUT:

        - ``key`` -- a key in ``self``

        OUTPUT:

        Either ``(None, 1)`` if the input ``key`` is not in the basis, or a tuple consisting of an UndirectedGraph and a sign.
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
        Return a string representation of ``self``.
        """
        return 'Basis consisting of representatives of isomorphism classes of labeled undirected graphs with no automorphisms that induce an odd permutation on edges'
