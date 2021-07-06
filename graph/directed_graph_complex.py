from .graph_complex import GraphCochain_dict, GraphComplex_dict, GraphCochain_vector, GraphComplex_vector
from .directed_graph_vector import DirectedGraphVector, DirectedGraphModule, DirectedGraphVector_dict, DirectedGraphModule_dict, DirectedGraphVector_vector, DirectedGraphModule_vector
from .directed_graph_basis import DirectedGraphComplexBasis
from util.misc import keydefaultdict
from functools import partial

class DirectedGraphCochain(DirectedGraphVector):
    """
    Cochain of a DirectedGraphComplex.
    """
    pass

class DirectedGraphComplex_(DirectedGraphModule):
    """
    Directed graph complex.
    """
    pass

class DirectedGraphCochain_dict(DirectedGraphCochain, DirectedGraphVector_dict, GraphCochain_dict):
    """
    Cochain of a DirectedGraphComplex (stored as a dictionary).
    """
    def __init__(self, parent, vector):
        """
        Initialize this graph cochain.
        """
        if not isinstance(parent, DirectedGraphComplex_dict):
            raise ValueError("parent must be a DirectedGraphComplex_dict")
        super().__init__(parent, vector)

class DirectedGraphComplex_dict(DirectedGraphComplex_, DirectedGraphModule_dict, GraphComplex_dict):
    """
    Directed graph complex (with elements stored as dictionaries).
    """
    def __init__(self, base_ring, connected=None, biconnected=None, min_degree=0, loops=True):
        """
        Initialize this graph complex.
        """
        graph_basis = DirectedGraphComplexBasis(connected=connected, biconnected=biconnected, min_degree=min_degree, loops=loops)
        super().__init__(base_ring, graph_basis)
        self.element_class = DirectedGraphCochain_dict

    def __repr__(self):
        """
        Return a string representation of this graph complex.
        """
        return 'Directed graph complex over {} with {}'.format(self._base_ring, self._graph_basis)

class DirectedGraphCochain_vector(DirectedGraphCochain, DirectedGraphVector_vector, GraphCochain_vector):
    """
    Cochain of a DirectedGraphComplex (stored as a dictionary of vectors).
    """
    def __init__(self, parent, vector):
        """
        Initialize this graph cochain.
        """
        if not isinstance(parent, DirectedGraphComplex_vector):
            raise ValueError("parent must be a DirectedGraphComplex_vector")
        super().__init__(parent, vector)

class DirectedGraphComplex_vector(DirectedGraphComplex_, DirectedGraphModule_vector, GraphComplex_vector):
    """
    Directed graph complex (with elements stored as dictionaries of vectors).
    """
    def __init__(self, base_ring, vector_constructor, matrix_constructor, connected=None, biconnected=None, min_degree=0, loops=True):
        """
        Initialize this graph complex.
        """
        if vector_constructor is None:
            raise ValueError('vector_constructor is required')
        if matrix_constructor is None:
            raise ValueError('matrix_constructor is required')
        graph_basis = DirectedGraphComplexBasis(connected=connected, biconnected=biconnected, min_degree=min_degree, loops=loops)
        super().__init__(base_ring, graph_basis, vector_constructor, matrix_constructor)
        self.element_class = DirectedGraphCochain_vector

    def __repr__(self):
        """
        Return a string representation of this graph complex.
        """
        return 'Directed graph complex over {} with {}'.format(self._base_ring, self._graph_basis)

def DirectedGraphComplex(base_ring, connected=None, biconnected=None, min_degree=0, loops=True, implementation='dict', vector_constructor=None, matrix_constructor=None):
    """
    Return the directed graph complex over ``base_ring`` with the given properties.
    """
    if implementation == 'dict':
        return DirectedGraphComplex_dict(base_ring, connected=connected, biconnected=biconnected, min_degree=min_degree, loops=loops)
    elif implementation == 'vector':
        return DirectedGraphComplex_vector(base_ring, vector_constructor, matrix_constructor, connected=connected, biconnected=biconnected, min_degree=min_degree, loops=loops)

