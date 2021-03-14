from .undirected_graph import UndirectedGraph
from .undirected_graph_vector_dict import UndirectedGraphVector_dict, UndirectedGraphModule_dict
from .undirected_graph_vector_vector import UndirectedGraphVector_vector, UndirectedGraphModule_vector
from .undirected_graph_basis import UndirectedGraphComplexBasis

class UndirectedGraphCochain_dict(UndirectedGraphVector_dict):
    """
    Cochain of an UndirectedGraphComplex (stored as a dictionary).
    """
    def __init__(self, parent, vector):
        """
        Initialize this graph cochain.
        """
        assert isinstance(parent, UndirectedGraphComplex_dict)
        super().__init__(parent, vector)

    def bracket(self, other):
        """
        Return the graph Lie bracket of this graph cochain with ``other``.
        """
        # TODO: optimize
        return sum(sum(self.homogeneous_part(v,e).insertion(i, other) for i in range(v)) for (v,e) in self.bi_gradings()) + sum((1 if e % 2 == 1 and f % 2 == 1 else -1)*sum(other.homogeneous_part(v,e).insertion(i, self.homogeneous_part(w,f)) for i in range(v)) for (v,e) in other.bi_gradings() for (w,f) in self.bi_gradings())

    def differential(self):
        """
        Return the graph differential of this graph cochain.
        """
        # TODO: optimize
        stick = self._parent(UndirectedGraph(2,[(0,1)]))
        return stick.bracket(self)

class UndirectedGraphComplex_dict(UndirectedGraphModule_dict):
    """
    Undirected graph complex (with elements stored as dictionaries).
    """
    def __init__(self, base_ring, connected=None, biconnected=None, min_degree=0):
        """
        Initialize this graph complex.
        """
        graph_basis = UndirectedGraphComplexBasis(connected=connected, biconnected=biconnected, min_degree=min_degree)
        super().__init__(base_ring, graph_basis)
        self.element_class = UndirectedGraphCochain_dict

    def __repr__(self):
        """
        Return a string representation of this graph complex.
        """
        return 'Undirected graph complex over {} with {}'.format(self._base_ring, self._graph_basis)

class UndirectedGraphCochain_vector(UndirectedGraphVector_vector):
    """
    Cochain of an UndirectedGraphComplex (stored as a dictionary of vectors).
    """
    def __init__(self, parent, vector):
        """
        Initialize this graph cochain.
        """
        assert isinstance(parent, UndirectedGraphComplex_vector)
        super().__init__(parent, vector)

    def bracket(self, other):
        """
        Return the graph Lie bracket of this graph cochain with ``other``.
        """
        # TODO: optimize
        return sum(sum(self.homogeneous_part(v,e).insertion(i, other) for i in range(v)) for (v,e) in self.bi_gradings()) + sum((1 if e % 2 == 1 and f % 2 == 1 else -1)*sum(other.homogeneous_part(v,e).insertion(i, self.homogeneous_part(w,f)) for i in range(v)) for (v,e) in other.bi_gradings() for (w,f) in self.bi_gradings())

    def differential(self):
        """
        Return the graph differential of this graph cochain.
        """
        # TODO: cache and optimize
        stick = self._parent(UndirectedGraph(2,[(0,1)]))
        return stick.bracket(self)

class UndirectedGraphComplex_vector(UndirectedGraphModule_vector):
    """
    Undirected graph complex (with elements stored as dictionaries of vectors).
    """
    def __init__(self, base_ring, vector_constructor, matrix_constructor, connected=None, biconnected=None, min_degree=0):
        """
        Initialize this graph complex.
        """
        graph_basis = UndirectedGraphComplexBasis(connected=connected, biconnected=biconnected, min_degree=min_degree)
        super().__init__(base_ring, graph_basis, vector_constructor, matrix_constructor)
        self.element_class = UndirectedGraphCochain_vector

    def __repr__(self):
        """
        Return a string representation of this graph complex.
        """
        return 'Undirected graph complex over {} with {}'.format(self._base_ring, self._graph_basis)

def UndirectedGraphComplex(base_ring, connected=None, biconnected=None, min_degree=0, implementation='dict', vector_constructor=None, matrix_constructor=None):
    """
    Return the undirected graph complex over ``base_ring`` with the given properties.
    """
    if implementation == 'dict':
        return UndirectedGraphComplex_dict(base_ring, connected=connected, biconnected=biconnected, min_degree=min_degree)
    elif implementation == 'vector':
        return UndirectedGraphComplex_vector(base_ring, vector_constructor, matrix_constructor, connected=connected, biconnected=biconnected, min_degree=min_degree)
