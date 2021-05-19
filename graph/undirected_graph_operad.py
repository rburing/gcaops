from .undirected_graph_vector import UndirectedGraphVector_dict, UndirectedGraphModule_dict
from .undirected_graph_basis import UndirectedGraphOperadBasis

class UndirectedGraphOperation_dict(UndirectedGraphVector_dict):
    """
    Element of an UndirectedGraphOperad (stored as a dictionary).
    """
    def __init__(self, parent, vector):
        """
        Initialize this graph operation.
        """
        if not isinstance(parent, UndirectedGraphOperad_dict):
            raise ValueError("parent must be a UndirectedGraphOperad_dict")
        super().__init__(parent, vector)

class UndirectedGraphOperad_dict(UndirectedGraphModule_dict):
    """
    Operad of undirected graphs (with elements stored as dictionaries).
    """
    def __init__(self, base_ring):
        """
        Initialize this graph operad.
        """
        graph_basis = UndirectedGraphOperadBasis()
        super().__init__(base_ring, graph_basis)
        self.element_class = UndirectedGraphOperation_dict

    def __repr__(self):
        """
        Return a string representation of this graph operad.
        """
        return 'Operad of undirected graphs over {}'.format(self._base_ring)

def UndirectedGraphOperad(base_ring):
    """
    Return the operad of undirected graphs over the given ``base_ring``.
    """
    return UndirectedGraphOperad_dict(base_ring)
