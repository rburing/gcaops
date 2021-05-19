from graph.undirected_graph import UndirectedGraph
from graph.graph_basis import GraphBasis
from util.misc import keydefaultdict
from util.undirected_graph_sage import undirected_graph_canonicalize, undirected_graph_generate, undirected_graph_has_odd_automorphism
from functools import partial

class UndirectedGraphBasis(GraphBasis):
    """
    Basis of a module spanned by undirected graphs.
    
    A basis consists of keys ``(v,e,index,...)`` where ``(v,e,index)`` identifies the isomorphism class of the graph.
    """
    graph_class = UndirectedGraph

class UndirectedGraphComplexBasis(UndirectedGraphBasis):
    """
    Basis consisting of representatives of isomorphism classes of undirected graphs with no automorphisms that induce an odd permutation on edges
    """
    def __init__(self, connected=None, biconnected=None, min_degree=0):
        """
        Initialize this basis.
        """
        if not min_degree in [0, 3]:
            raise ValueError('min_degree can only be 0 or 3')
        self._connected = connected
        self._biconnected = biconnected
        self._min_degree = min_degree
        self._graphs = keydefaultdict(partial(__class__._generate_graphs, self))

    def _generate_graphs(self, bi_grading):
        """
        Return a list of all the graphs in this basis in the given ``bi_grading``.
        """
        v, e = bi_grading
        graphs = []
        for g in undirected_graph_generate(v, e, connected=self._connected, biconnected=self._biconnected, min_degree=self._min_degree):
            if not undirected_graph_has_odd_automorphism(g):
                graphs.append(g)
        return graphs

    def graph_to_key(self, graph):
        """
        Return a tuple consisting of the key in this basis and the sign factor such that ``graph`` equals the sign times the graph identified by the key.

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

        - ``key`` -- a key in this basis

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
        Return a string representation of this basis.
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

class UndirectedGraphOperadBasis(UndirectedGraphBasis):
    """
    Basis consisting of representatives of isomorphism classes of labeled undirected graphs with no automorphisms that induce an odd permutation on edges
    """
    def __init__(self):
        """
        Initialize this basis.
        """
        self._graphs = keydefaultdict(partial(__class__._generate_graphs, self))

    def _generate_graphs(self, bi_grading):
        """
        Return a list of all the graphs in this basis in the given ``bi_grading``.
        """
        v, e = bi_grading
        graphs = []
        for g in undirected_graph_generate(v, e):
            if not undirected_graph_has_odd_automorphism(g):
                graphs.append(g)
        return graphs

    def graph_to_key(self, graph):
        """
        Return a tuple consisting of the key in this basis and the sign factor such that ``graph`` equals the sign times the graph identified by the key.

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

        - ``key`` -- a key in this basis

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
        Return a string representation of this basis.
        """
        return 'Basis consisting of representatives of isomorphism classes of labeled undirected graphs with no automorphisms that induce an odd permutation on edges'

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
