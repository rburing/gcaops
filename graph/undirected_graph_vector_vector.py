from .undirected_graph import UndirectedGraph
from .undirected_graph_basis import UndirectedGraphBasis
from util.misc import keydefaultdict
from itertools import product

class UndirectedGraphVector_vector:
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
        assert isinstance(parent, UndirectedGraphModule_vector)
        self._parent = parent
        self._vectors = keydefaultdict(lambda key: self._parent._vector_constructor(self._parent.base_ring(), self._parent.basis().cardinality(*key)))
        for bi_grading in vectors:
            self._vectors[bi_grading] = vectors[bi_grading]

    def __repr__(self):
        """
        Return a string representation of this graph vector.
        """
        terms = []
        for (bi_grading, vector) in self._vectors.items():
            for (k,c) in vector.items():
                g, sign = self._parent._graph_basis.key_to_graph(bi_grading + (k,))
                c *= sign
                c_str = repr(c)
                if c_str != '1':
                    c_str = '({})'.format(c)
                terms.append('{}*{}'.format(c_str, repr(g)))
        if len(terms) > 0:
            return ' + '.join(terms)
        else:
            return '0'

    def parent(self):
        """
        Return the parent UndirectedGraphModule that this graph vector belongs to.
        """
        return self._parent

    def copy(self):
        """
        Return a copy of this graph vector.
        """
        return self.__class__(self._parent, self._vectors)

    __pos__ = copy

    def __neg__(self):
        """
        Return the negative of this graph vector.
        """
        return self.__class__(self._parent, {k : -v for (k,v) in self._vectors.items()})

    def __add__(self, other):
        """
        Return this graph vector added to ``other``.
        """
        if isinstance(other, self.__class__):
            v = self._vectors.copy()
            for k in other._vectors:
                v[k] += other._vectors[k]
            return self.__class__(self._parent, v)
        elif other == 0:
            return self.copy()

    def __radd__(self, other):
        """
        Return ``other`` added to this graph vector.
        """
        return self + other

    def __sub__(self, other):
        """
        Return ``other`` subtracted from this graph vector.
        """
        if isinstance(other, self.__class__):
            v = self._vectors.copy()
            for k in other._vectors:
                v[k] -= other._vectors[k]
            return self.__class__(self._parent, v)
        elif other == 0:
            return self.copy()

    def __rsub__(self, other):
        """
        Return this graph vector subtracted from ``other``.
        """
        return -(self - other)

    def __mul__(self, other):
        """
        Return this graph vector multiplied by ``other``.
        """
        v = self._vectors.copy()
        if other in self._parent.base_ring():
            for k in v:
                v[k] *= other
            return self.__class__(self._parent, v)
        else:
            raise NotImplementedError

    def __rmul__(self, other):
        """
        Return ``other`` multiplied by this graph vector.
        """
        return self * other

    def __eq__(self, other):
        """
        Return ``True`` if this graph vector is equal to ``other`` and ``False`` otherwise.
        """
        difference = self - other
        for k in difference._vectors:
            if not difference._vectors[k].is_zero():
                return False
        return True

    def bi_gradings(self):
        """
        Return the set of tuples ``(v,e)`` such that this graph vector contains terms with ``v`` vertices and ``e`` edges.
        """
        return set(self._vectors.keys())

    def nvertices(self):
        """
        Return the number of vertices in each graph in this graph vector.

        ASSUMPTIONS:

        Assumes all graphs in this graph vector have the same number of vertices.
        """
        for bi_grading in self._vectors:
            if not self._vectors[bi_grading].is_zero():
                return bi_grading[0]

    def homogeneous_part(self, vertices, edges):
        """
        Return the homogeneous part of this graph vector consisting only of terms with the given number of ``vertices`` and ``edges``.
        """
        bi_grading = (vertices, edges)
        return self.__class__(self._parent, { bi_grading : self._vectors[bi_grading]})

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
                            term = UndirectedGraph(len(user) + len(victim) - 1, [tuple(e) for e in user_edges] + victim_edges)
                            terms.append([user_coeff*victim_coeff, term])
        return self._parent(terms)

class UndirectedGraphModule_vector:
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
        self._base_ring = base_ring
        if not isinstance(graph_basis, UndirectedGraphBasis):
            raise ValueError('graph_basis must be UndirectedGraphBasis')
        self._graph_basis = graph_basis
        self._vector_constructor = vector_constructor
        self._matrix_constructor = matrix_constructor
        self.element_class = UndirectedGraphVector_vector

    def base_ring(self):
        """
        Return the base ring of this module.
        """
        return self._base_ring

    def basis(self):
        """
        Return the basis of this module.
        """
        return self._graph_basis

    def __repr__(self):
        """
        Return a string representation of this module.
        """
        return 'Module over {} with {}'.format(self._base_ring, self._graph_basis)

    def zero(self):
        """
        Return the zero vector in this module.
        """
        return self.element_class(self, {})

    def __call__(self, arg):
        """
        Convert ``arg`` into an element of this module.
        """
        if isinstance(arg, UndirectedGraph):
            key, sign = self._graph_basis.graph_to_key(arg)
            if key is not None:
                bi_grading = key[:2]
                index = key[2]
                v = self.zero()
                v._vectors[bi_grading][index] = self.base_ring().one() * sign
                return v
            else: # must be zero
                return self.zero()
        elif isinstance(arg, list):
            v = self.element_class(self, {})
            for term in arg:
                coeff, graph = term
                key, sign = self._graph_basis.graph_to_key(graph)
                if key is not None:
                    coeff *= sign
                    bi_grading = key[:2]
                    index = key[2]
                    v._vectors[bi_grading][index] += coeff
            return v
        elif arg == 0:
            return self.zero()
