from .directed_graph_vector import DirectedGraphVector_dict, DirectedGraphModule_dict
from .directed_graph_basis import DirectedGraphOperadBasis

class DirectedGraphOperation_dict(DirectedGraphVector_dict):
    """
    Element of a DirectedGraphOperad (stored as a dictionary).
    """
    def __init__(self, parent, vector):
        """
        Initialize this graph operation.
        """
        if not isinstance(parent, DirectedGraphOperad_dict):
            raise ValueError("parent must be a DirectedGraphOperad_dict")
        super().__init__(parent, vector)

class DirectedGraphOperad_dict(DirectedGraphModule_dict):
    """
    Operad of directed graphs (with elements stored as dictionaries).
    """
    def __init__(self, base_ring):
        """
        Initialize this graph operad.
        """
        graph_basis = DirectedGraphOperadBasis()
        super().__init__(base_ring, graph_basis)
        self.element_class = DirectedGraphOperation_dict

    def __repr__(self):
        """
        Return a string representation of this graph operad.
        """
        return 'Operad of directed graphs over {}'.format(self._base_ring)

def DirectedGraphOperad(base_ring):
    """
    Return the operad of directed graphs over the given ``base_ring``.
    """
    return DirectedGraphOperad_dict(base_ring)

