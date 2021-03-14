from .undirected_graph import UndirectedGraph
from .undirected_graph_vector import UndirectedGraphVector, UndirectedGraphModule
from .undirected_graph_basis import UndirectedGraphComplexBasis

class UndirectedGraphCochain(UndirectedGraphVector):
    """
    Cochain of an UndirectedGraphComplex.
    """
    def __init__(self, parent, vector):
        """
        Initialize this graph cochain.
        """
        assert isinstance(parent, UndirectedGraphComplex)
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

class UndirectedGraphComplex(UndirectedGraphModule):
    """
    Undirected graph complex.
    """
    def __init__(self, base_ring, connected=None, biconnected=None, min_degree=0):
        """
        Initialize this graph complex.
        """
        graph_basis = UndirectedGraphComplexBasis(connected=connected, biconnected=biconnected, min_degree=min_degree)
        super().__init__(base_ring, graph_basis)
        self.element_class = UndirectedGraphCochain

    def __repr__(self):
        """
        Return a string representation of this graph complex.
        """
        return 'Undirected graph complex over {} with {}'.format(self._base_ring, self._graph_basis)
