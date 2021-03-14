from .undirected_graph_vector import UndirectedGraphVector, UndirectedGraphModule
from .undirected_graph_basis import UndirectedGraphOperadBasis

class UndirectedGraphOperation(UndirectedGraphVector):
    """
    Element of an UndirectedGraphOperad.
    """
    def __init__(self, parent, vector):
        """
        Initialize this graph operation.
        """
        assert isinstance(parent, UndirectedGraphOperad)
        super().__init__(parent, vector)

class UndirectedGraphOperad(UndirectedGraphModule):
    """
    Operad of undirected graphs.
    """
    def __init__(self, base_ring):
        """
        Initialize this graph operad.
        """
        graph_basis = UndirectedGraphOperadBasis()
        super().__init__(base_ring, graph_basis)
        self.element_class = UndirectedGraphOperation

    def __repr__(self):
        """
        Return a string representation of this graph operad.
        """
        return 'Operad of undirected graphs over {}'.format(self._base_ring)
