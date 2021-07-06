from .graph_complex import GraphCochain_dict, GraphComplex_dict, GraphCochain_vector, GraphComplex_vector
from .undirected_graph_vector import UndirectedGraphVector, UndirectedGraphModule, UndirectedGraphVector_dict, UndirectedGraphModule_dict, UndirectedGraphVector_vector, UndirectedGraphModule_vector
from .undirected_graph_basis import UndirectedGraphComplexBasis
from util.misc import keydefaultdict
from functools import partial
from abc import abstractmethod

class UndirectedGraphCochain(UndirectedGraphVector):
    """
    Cochain of an UndirectedGraphComplex.
    """
    @abstractmethod
    def _indices_and_coefficients(self, bi_grading):
        """
        Return an iterator over tuples ``(index, coefficient)`` in this graph cochain.
        """
        pass

class UndirectedGraphComplex_(UndirectedGraphModule):
    """
    Undirected graph complex.
    """
    pass

class UndirectedGraphCochain_dict(UndirectedGraphCochain, UndirectedGraphVector_dict, GraphCochain_dict):
    """
    Cochain of an UndirectedGraphComplex (stored as a dictionary).
    """
    def __init__(self, parent, vector):
        """
        Initialize this graph cochain.
        """
        if not isinstance(parent, UndirectedGraphComplex_dict):
            raise ValueError("parent must be a UndirectedGraphComplex_dict")
        super().__init__(parent, vector)

    def _indices_and_coefficients(self, bi_grading):
        """
        Return an iterator over tuples ``(index, coefficient)`` in this graph cochain.
        """
        for key in self._vector:
            c = self._vector[key]
            if c.is_zero():
                continue
            num_vertices, num_edges, index = key
            if (num_vertices, num_edges) == bi_grading:
                yield (index, c)

class UndirectedGraphComplex_dict(UndirectedGraphComplex_, UndirectedGraphModule_dict, GraphComplex_dict):
    """
    Undirected graph complex (with elements stored as dictionaries).
    """
    def __init__(self, base_ring, connected=None, biconnected=None, min_degree=0):
        """
        Initialize this graph complex.
        """
        if not min_degree in [0, 3]:
            raise ValueError('min_degree can only be 0 or 3')
        graph_basis = UndirectedGraphComplexBasis(connected=connected, biconnected=biconnected, min_degree=min_degree)
        super().__init__(base_ring, graph_basis)
        self.element_class = UndirectedGraphCochain_dict

    def __repr__(self):
        """
        Return a string representation of this graph complex.
        """
        return 'Undirected graph complex over {} with {}'.format(self._base_ring, self._graph_basis)

class UndirectedGraphCochain_vector(UndirectedGraphCochain, UndirectedGraphVector_vector, GraphCochain_vector):
    """
    Cochain of an UndirectedGraphComplex (stored as a dictionary of vectors).
    """
    def __init__(self, parent, vector):
        """
        Initialize this graph cochain.
        """
        if not isinstance(parent, UndirectedGraphComplex_vector):
            raise ValueError("parent must be a UndirectedGraphComplex_vector")
        super().__init__(parent, vector)

    def _indices_and_coefficients(self, bi_grading):
        """
        Return an iterator over tuples ``(index, coefficient)`` in this graph cochain.
        """
        yield from self._vectors[bi_grading].items()

class UndirectedGraphComplex_vector(UndirectedGraphComplex_, UndirectedGraphModule_vector, GraphComplex_vector):
    """
    Undirected graph complex (with elements stored as dictionaries of vectors).
    """
    def __init__(self, base_ring, vector_constructor, matrix_constructor, connected=None, biconnected=None, min_degree=0):
        """
        Initialize this graph complex.
        """
        if vector_constructor is None:
            raise ValueError('vector_constructor is required')
        if matrix_constructor is None:
            raise ValueError('matrix_constructor is required')
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
