from abc import abstractmethod
from copy import copy
from .graph_vector import GraphVector, GraphModule
from .graph_vector_dict import GraphVector_dict, GraphModule_dict
from .graph_vector_vector import GraphVector_vector, GraphModule_vector
from .directed_graph_basis import DirectedGraphBasis
# for conversion:
from .undirected_graph import UndirectedGraph
from .undirected_graph_vector import UndirectedGraphVector

class DirectedGraphVector(GraphVector):
    """
    Vector representing a linear combination of directed graphs.
    """
    @abstractmethod
    def filter(self, max_out_degree=None):
        """
        Return the graph vector which is the summand of this graph vector containing exactly those graphs that pass the filter.
        """
        pass

class DirectedGraphModule(GraphModule):
    """
    Module spanned by directed graphs.
    """
    def __call__(self, arg):
        """
        Convert ``arg`` into an element of this module.
        """
        # NOTE: the calls to super().__call__ here will go to the concrete implementation
        if isinstance(arg, UndirectedGraph):
            result = self.zero()
            for g in arg.orientations():
                result += super().__call__(g)
            return result
        elif isinstance(arg, UndirectedGraphVector):
            result = self.zero()
            for (c,g) in arg:
                for h in g.orientations():
                    result += c*super().__call__(h)
            return result
        else:
            return super().__call__(arg)

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

    def filter(self, max_out_degree=None):
        """
        Return the graph vector which is the summand of this graph vector containing exactly those graphs that pass the filter.
        """
        new_vector = {}
        for key in self._vector:
            c = self._vector[key]
            if c.is_zero():
                continue
            g, sign = self._parent._graph_basis.key_to_graph(key)
            if any(d > max_out_degree for d in g.out_degrees()):
                c = self._parent.base_ring().zero()
            new_vector[key] = c
        return self.__class__(self._parent, new_vector)

class DirectedGraphModule_dict(DirectedGraphModule, GraphModule_dict):
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

    def filter(self, max_out_degree=None):
        """
        Return the graph vector which is the summand of this graph vector containing exactly those graphs that pass the filter.
        """
        v = {}
        for (bi_grading, vector) in self._vectors.items():
            v[bi_grading] = copy(vector)
            for j in vector.nonzero_positions():
                g, sign = self._parent._graph_basis.key_to_graph(bi_grading + (j,))
                if any(d > max_out_degree for d in g.out_degrees()):
                    v[bi_grading][j] = self._parent.base_ring().zero()
        return self.__class__(self._parent, v)

class DirectedGraphModule_vector(DirectedGraphModule, GraphModule_vector):
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
