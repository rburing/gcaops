from .graph_vector_vector import GraphVector_vector, GraphModule_vector
from .undirected_graph_basis import UndirectedGraphBasis

class UndirectedGraphVector_vector(GraphVector_vector):
    """
    Vector representing a linear combination of undirected graphs (stored as a dictionary of vectors).
    """
    def __init__(self, parent, vectors):
        """
        Initialize this graph vector.

        INPUT:

        - ``parent`` -- an UndirectedGraphModule

        - ``vectors`` -- a dictionary, mapping bi-gradings to sparse vectors of coefficients with respect to the basis of ``parent``
        """
        if not isinstance(parent, UndirectedGraphModule_vector):
            raise ValueError("parent must be a UndirectedGraphModule_vector")
        super().__init__(parent, vectors)

class UndirectedGraphModule_vector(GraphModule_vector):
    """
    Module spanned by undirected graphs (with elements stored as dictionaries of vectors).
    """
    def __init__(self, base_ring, graph_basis, vector_constructor, matrix_constructor):
        """
        Initialize this undirected graph module.

        INPUT:

        - ``base_ring`` -- a ring, to be used as the ring of coefficients

        - ``graph_basis`` -- an UndirectedGraphBasis

        - ``vector_constructor`` -- constructor of (sparse) vectors

        - ``matrix_constructor`` -- constructor of (sparse) matrices
        """
        if not isinstance(graph_basis, UndirectedGraphBasis):
            raise ValueError('graph_basis must be an UndirectedGraphBasis')
        super().__init__(base_ring, graph_basis, vector_constructor, matrix_constructor)
        self.element_class = UndirectedGraphVector_vector
