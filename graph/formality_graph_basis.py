from graph.formality_graph import FormalityGraph
from graph.graph_basis import GraphBasis
from graph.graph_cache import formality_graph_cache
from util.misc import keydefaultdict
from functools import partial

class FormalityGraphBasis(GraphBasis):
    """
    Basis of a module spanned by formality graphs.

    A basis consists of keys ``(gv,av,e,index,...)`` where ``(gv,av,e,index)`` identifies the isomorphism class of the graph.
    """
    graph_class = FormalityGraph
    grading_size = 3

class FormalityGraphComplexBasis(FormalityGraphBasis):
    """
    Basis consisting of representatives of isomorphism classes of formality graphs with no automorphisms that induce an odd permutation on edges.
    """
    def __init__(self, positive_differential_order=None):
        """
        Initialize this basis.
        """
        self._positive_differential_order = positive_differential_order
        self._graphs = keydefaultdict(partial(formality_graph_cache.graphs, positive_differential_order=positive_differential_order, has_odd_automorphism=False))

    def graph_to_key(self, graph):
        """
        Return a tuple consisting of the key in this basis and the sign factor such that ``graph`` equals the sign times the graph identified by the key.

        INPUT:

        - ``graph`` -- a FormalityGraph

        OUTPUT:

        Either ``(None, 1)`` if the input ``graph`` is not in the span of the basis, or a tuple consisting of a key and a sign, where a key is a tuple consisting of the number of ground vertices, the number of aerial vertices, the number of edges, and the index of the graph in the list.
        """
        g, _, sign = formality_graph_cache.canonicalize_graph(graph)
        gv, av, e = g.num_ground_vertices(), g.num_aerial_vertices(), len(g.edges())
        try:
            index = self._graphs[gv,av,e].index(g)
            return (gv,av,e,index), sign
        except ValueError:
            return None, 1

    def key_to_graph(self, key):
        """
        Return a tuple consisting of a FormalityGraph and the sign factor such that the sign times the graph equals the graph identified by the key.

        INPUT:

        - ``key`` -- a key in this basis

        OUTPUT:

        Either ``(None, 1)`` if the input ``key`` is not in the basis, or a tuple consisting of a FormalityGraph and a sign which is always +1.
        """
        gv, av, e, index = key 
        try:
            return self._graphs[gv,av,e][index], 1
        except IndexError:
            return None, 1

    def __repr__(self):
        """
        Return a string representation of this basis.
        """
        filters = []
        if self._positive_differential_order:
            filters.append('of positive differential order')
        if filters:
            filters_str = ' ({})'.format(', '.join(filters))
        else:
            filters_str = ''
        return 'Basis consisting of representatives of isomorphism classes of formality graphs{} with no automorphisms that induce an odd permutation on edges'.format(filters_str)

    def graph_properties(self):
        """
        Return a dictionary containing the properties of the graphs in this basis.
        """
        return {'positive_differential_order' : self._positive_differential_order, 'has_odd_automorphism' : False}

    def graphs(self, num_ground_vertices, num_aerial_vertices, num_edges):
        """
        Return the list of graphs in this basis with the given ``num_ground_vertices``, ``num_aerial_vertices`` and ``num_edges``.
        """
        return self._graphs[num_ground_vertices, num_aerial_vertices, num_edges]

    def cardinality(self, num_ground_vertices, num_aerial_vertices, num_edges):
        """
        Return the number of graphs in this basis with the given ``num_ground_vertices``, ``num_aerial_vertices`` and ``num_edges``.
        """
        return len(self._graphs[num_ground_vertices, num_aerial_vertices, num_edges])

class FormalityGraphOperadBasis(FormalityGraphBasis):
    """
    Basis consisting of labeled formality graphs with no automorphisms that induce an odd permutation on edges
    """
    def __init__(self, positive_differential_order=None):
        """
        Initialize this basis.
        """
        self._positive_differential_order = positive_differential_order
        self._graphs = keydefaultdict(partial(formality_graph_cache.graphs, positive_differential_order=positive_differential_order, has_odd_automorphism=False))

    def graph_to_key(self, graph):
        """
        Return a tuple consisting of the key in this basis and the sign factor such that ``graph`` equals the sign times the graph identified by the key.

        INPUT:

        - ``graph`` -- a FormalityGraph

        OUTPUT:

        Either ``(None, 1)`` if the input ``graph`` is not in the span of the basis, or a tuple consisting of a key and a sign, where a key is a tuple consisting of the number of ground vertices, the number of aerial vertices, the number of edges, the index of the graph in the list, followed by a permutation of vertices.
        """
        g, undo_canonicalize, sign = formality_graph_cache.canonicalize_graph(graph)
        gv, av, e = g.num_ground_vertices(), g.num_aerial_vertices(), len(g.edges())
        try:
            index = self._graphs[gv,av,e].index(g)
            return (gv,av,e,index) + tuple(undo_canonicalize), sign
        except ValueError:
            return None, 1

    def key_to_graph(self, key):
        """
        Return a tuple consisting of a FormalityGraph and the sign factor such that the sign times the graph equals the graph identified by the key.

        INPUT:

        - ``key`` -- a key in this basis

        OUTPUT:

        Either ``(None, 1)`` if the input ``key`` is not in the basis, or a tuple consisting of a FormalityGraph and a sign.
        """
        gv, av, e, index = key[:4]
        undo_canonicalize = key[4:]
        try:
            G = self._graphs[gv,av,e][index]
            g = G.relabeled(undo_canonicalize)
            sign = g.canonicalize_edges()
            return g, sign
        except IndexError:
            return None, 1

    def __repr__(self):
        """
        Return a string representation of this basis.
        """
        filters = []
        if self._positive_differential_order:
            filters.append('of positive differential order')
        if filters:
            filters_str = ' ({})'.format(', '.join(filters))
        else:
            filters_str = ''
        return 'Basis consisting of labeled formality graphs{} with no automorphisms that induce an odd permutation on edges'.format(filters_str)

    def graph_properties(self):
        """
        Return a dictionary containing the properties of the graphs in this basis.
        """
        return {'positive_differential_order' : self._positive_differential_order, 'has_odd_automorphism' : False}

