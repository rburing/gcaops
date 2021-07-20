from itertools import product
from .graph_vector import GraphVector, GraphModule
from .graph_vector_dict import GraphVector_dict, GraphModule_dict
from .graph_vector_vector import GraphVector_vector, GraphModule_vector
from .undirected_graph_basis import UndirectedGraphBasis

class UndirectedGraphVector(GraphVector):
    """
    Vector representing a linear combination of undirected graphs.
    """
    pass

class UndirectedGraphModule(GraphModule):
    """
    Module spanned by undirected graphs.
    """
    pass

class UndirectedGraphVector_dict(UndirectedGraphVector, GraphVector_dict):
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

    def nvertices(self):
        """
        Return the number of vertices in each graph in this graph vector.

        ASSUMPTIONS:

        Assumes all graphs in this graph vector have the same number of vertices.
        """
        for key in self._vector:
            v, e = key[:2]
            if not self._vector[key].is_zero():
                return v

    def nedges(self):
        """
        Return the number of edges in each graph in this graph vector.

        ASSUMPTIONS:

        Assumes all graphs in this graph vector have the same number of edges.
        """
        for key in self._vector:
            v, e = key[:2]
            if not self._vector[key].is_zero():
                return e

    def insertion(self, position, other):
        """
        Return the insertion of ``other`` into this graph vector at the vertex ``position``.
        """
        # TODO: cache when self and other are in normal form. when not, use symmetric group action + operad axioms to deduce result.
        terms = []
        for user_key in self._vector:
            user_coeff = self._vector[user_key]
            if user_coeff.is_zero():
                continue
            for victim_key in other._vector:
                victim_coeff = other._vector[victim_key]
                if victim_coeff.is_zero():
                    continue
                user, user_sign = self._parent._graph_basis.key_to_graph(user_key)
                user_coeff *= user_sign
                victim, victim_sign = other._parent._graph_basis.key_to_graph(victim_key)
                victim_coeff *= victim_sign
                # relabel user (vertices > position are shifted to make room for victim)
                user_edges = [[a + len(victim) - 1 if a > position else a, b + len(victim) - 1 if b > position else b] for (a,b) in user.edges()]
                # relabel victim
                victim_edges = [(position + a, position + b) for (a,b) in victim.edges()]
                # find edges which are incident to position
                incident = [(i,user_edges[i].index(position)) for i in range(len(user_edges)) if position in user_edges[i]]
                # loop over all possible new endpoints (in victim) for these edges
                for endpoints in product(range(len(victim)), repeat=len(incident)):
                    # redirect edges (which were incident to position) to victim
                    for k in range(len(incident)):
                        a, b = incident[k]
                        user_edges[a][b] = position + endpoints[k]
                    # NOTE: the convention is that victim edges go last:
                    term = self._parent._graph_basis.graph_class(len(user) + len(victim) - 1, [tuple(e) for e in user_edges] + victim_edges)
                    terms.append([user_coeff*victim_coeff, term])
        return self._parent(terms)

class UndirectedGraphModule_dict(UndirectedGraphModule, GraphModule_dict):
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

class UndirectedGraphVector_vector(UndirectedGraphVector, GraphVector_vector):
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

    def nvertices(self):
        """
        Return the number of vertices in each graph in this graph vector.

        ASSUMPTIONS:

        Assumes all graphs in this graph vector have the same number of vertices.
        """
        for bi_grading in self._vectors:
            if not self._vectors[bi_grading].is_zero():
                return bi_grading[0]

    def nedges(self):
        """
        Return the number of edges in each graph in this graph vector.

        ASSUMPTIONS:

        Assumes all graphs in this graph vector have the same number of edges.
        """
        for bi_grading in self._vectors:
            if not self._vectors[bi_grading].is_zero():
                return bi_grading[1]

    def insertion(self, position, other):
        """
        Return the insertion of ``other`` into this graph vector at the vertex ``position``.
        """
        terms = []
        for (user_bigrading, user_vector) in self._vectors.items():
            for (user_idx, user_coeff) in user_vector.items():
                user_key = user_bigrading + (user_idx,)
                user, user_sign = self._parent._graph_basis.key_to_graph(user_key)
                user_coeff *= user_sign
                for (victim_bigrading, victim_vector) in other._vectors.items():
                    for (victim_idx, victim_coeff) in victim_vector.items():
                        victim_key = victim_bigrading + (victim_idx,)
                        victim, victim_sign = other._parent._graph_basis.key_to_graph(victim_key)
                        victim_coeff *= victim_sign
                        # relabel user (vertices > position are shifted to make room for victim)
                        user_edges = [[a + len(victim) - 1 if a > position else a, b + len(victim) - 1 if b > position else b] for (a,b) in user.edges()]
                        # relabel victim
                        victim_edges = [(position + a, position + b) for (a,b) in victim.edges()]
                        # find edges which are incident to position
                        incident = [(i,user_edges[i].index(position)) for i in range(len(user_edges)) if position in user_edges[i]]
                        # loop over all possible new endpoints (in victim) for these edges
                        for endpoints in product(range(len(victim)), repeat=len(incident)):
                            # redirect edges (which were incident to position) to victim
                            for k in range(len(incident)):
                                a, b = incident[k]
                                user_edges[a][b] = position + endpoints[k]
                            # NOTE: the convention is that victim edges go last:
                            term = self._parent._graph_basis.graph_class(len(user) + len(victim) - 1, [tuple(e) for e in user_edges] + victim_edges)
                            terms.append([user_coeff*victim_coeff, term])
        return self._parent(terms)

class UndirectedGraphModule_vector(UndirectedGraphModule, GraphModule_vector):
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
