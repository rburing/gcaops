from .graph_vector_dict import GraphVector_dict, GraphModule_dict
from .undirected_graph_basis import UndirectedGraphBasis

class UndirectedGraphVector_dict(GraphVector_dict):
    """
    Vector representing a linear combination of undirected graphs (stored as a dictionary).
    """
    def __init__(self, parent, vector):
        """
        Initialize this undirected graph vector.

        INPUT:

        - ``parent`` -- an UndirectedGraphModule

        - ``vector`` -- a dictionary, representing a sparse vector of coefficients with respect to the basis of ``parent``
        """
        if not isinstance(parent, UndirectedGraphModule_dict):
            raise ValueError("parent must be a UndirectedGraphModule_dict")
        super().__init__(parent, vector)

class UndirectedGraphModule_dict(GraphModule_dict):
    """
    Module spanned by undirected graphs (with elements stored as dictionaries).
    """
    def __init__(self, base_ring, graph_basis):
        """
        Initialize this undirected graph module.

        INPUT:

        - ``base_ring`` -- a ring, to be used as the ring of coefficients

        - ``graph_basis`` -- an UndirectedGraphBasis
        """
        if not isinstance(graph_basis, UndirectedGraphBasis):
            raise ValueError('graph_basis must be an UndirectedGraphBasis')
        super().__init__(base_ring, graph_basis)
        self.element_class = UndirectedGraphVector_dict
