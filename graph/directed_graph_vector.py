from .graph_vector import GraphVector
from .graph_vector_dict import GraphVector_dict, GraphModule_dict
from .graph_vector_vector import GraphVector_vector, GraphModule_vector
from .directed_graph_basis import DirectedGraphBasis

class DirectedGraphVector(GraphVector):
    """
    Vector representing a linear combination of directed graphs.
    """
    pass

class DirectedGraphVector_dict(DirectedGraphVector, GraphVector_dict):
    """
    Vector representing a linear combination of directed graphs (stored as a dictionary).
    """
    def __init__(self, parent, vector):
        """
        Initialize this directed graph vector.

        INPUT:

        - ``parent`` -- a DirectedGraphModule

        - ``vector`` -- a dictionary, representing a sparse vector of coefficients with respect to the basis of ``parent``
        """
        if not isinstance(parent, DirectedGraphModule_dict):
            raise ValueError("parent must be a DirectedGraphModule_dict")
        super().__init__(parent, vector)

class DirectedGraphModule_dict(GraphModule_dict):
    """
    Module spanned by directed graphs (with elements stored as dictionaries).
    """
    def __init__(self, base_ring, graph_basis):
        """
        Initialize this directed graph module.

        INPUT:

        - ``base_ring`` -- a ring, to be used as the ring of coefficients

        - ``graph_basis`` -- a DirectedGraphBasis
        """
        if not isinstance(graph_basis, DirectedGraphBasis):
            raise ValueError('graph_basis must be a DirectedGraphBasis')
        super().__init__(base_ring, graph_basis)
        self.element_class = DirectedGraphVector_dict

class DirectedGraphVector_vector(DirectedGraphVector, GraphVector_vector):
    """
    Vector representing a linear combination of directed graphs (stored as a dictionary of vectors).
    """
    def __init__(self, parent, vectors):
        """
        Initialize this graph vector.

        INPUT:

        - ``parent`` -- a DirectedGraphModule

        - ``vectors`` -- a dictionary, mapping bi-gradings to sparse vectors of coefficients with respect to the basis of ``parent``
        """
        if not isinstance(parent, DirectedGraphModule_vector):
            raise ValueError("parent must be a DirectedGraphModule_vector")
        super().__init__(parent, vectors)

class DirectedGraphModule_vector(GraphModule_vector):
    """
    Module spanned by directed graphs (with elements stored as dictionaries of vectors).
    """
    def __init__(self, base_ring, graph_basis, vector_constructor, matrix_constructor):
        """
        Initialize this directed graph module.

        INPUT:

        - ``base_ring`` -- a ring, to be used as the ring of coefficients

        - ``graph_basis`` -- a DirectedGraphBasis

        - ``vector_constructor`` -- constructor of (sparse) vectors

        - ``matrix_constructor`` -- constructor of (sparse) matrices
        """
        if not isinstance(graph_basis, DirectedGraphBasis):
            raise ValueError('graph_basis must be a DirectedGraphBasis')
        super().__init__(base_ring, graph_basis, vector_constructor, matrix_constructor)
        self.element_class = DirectedGraphVector_vector
