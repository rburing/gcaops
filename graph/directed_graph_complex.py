from .graph_complex import GraphCochain_dict, GraphComplex_dict, GraphCochain_vector, GraphComplex_vector
from .directed_graph_vector import DirectedGraphVector, DirectedGraphModule, DirectedGraphVector_dict, DirectedGraphModule_dict, DirectedGraphVector_vector, DirectedGraphModule_vector
from .directed_graph_basis import DirectedGraphComplexBasis
from .undirected_graph_complex import UndirectedGraphCochain
from util.misc import keydefaultdict
from functools import partial
from abc import abstractmethod

class DirectedGraphCochain(DirectedGraphVector):
    """
    Cochain of a DirectedGraphComplex.
    """
    @abstractmethod
    def _add_to_coeff_by_index(self, bi_grading, index, summand):
        """
        Add ``summand`` to the specified coefficient in this graph cochain.
        """
        pass

class DirectedGraphComplex_(DirectedGraphModule):
    """
    Directed graph complex.
    """
    def __call__(self, arg):
        if isinstance(arg, UndirectedGraphCochain):
            # We try to use the cache, if graph properties are compatible
            undirected_properties = arg.parent().basis().graph_properties()
            directed_properties = self.basis().graph_properties()
            del directed_properties['loops']
            if directed_properties != undirected_properties:
                return super().__call__(arg) # NOTE: falling back to slow implementation
            result = self.zero()
            for bi_grading in arg.bi_gradings():
                for (undirected_graph_idx, coefficient) in arg._indices_and_coefficients(bi_grading):
                    for (directed_graph_idx, coefficient2) in self._graph_basis._undirected_to_directed_coeffs(bi_grading, undirected_graph_idx):
                        result._add_to_coeff_by_index(bi_grading, directed_graph_idx, coefficient * coefficient2)
            return result
        else:
            return super().__call__(arg)

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

    def _add_to_coeff_by_index(self, bi_grading, index, summand):
        """
        Add ``summand`` to the specified coefficient in this graph cochain.
        """
        self._vector[bi_grading + (index,)] += summand

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

    def _add_to_coeff_by_index(self, bi_grading, index, summand):
        """
        Add ``summand`` to the specified coefficient in this graph cochain.
        """
        self._vectors[bi_grading][index] += summand

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

