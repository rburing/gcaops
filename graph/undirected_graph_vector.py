from .undirected_graph import UndirectedGraph
from .undirected_graph_basis import UndirectedGraphBasis
from collections import defaultdict
from itertools import product

class UndirectedGraphVector:
    """
    Vector representing a linear combination of undirected graphs.
    """
    def __init__(self, parent, vector):
        """
        Initialize ``self``.

        INPUT:

        - ``parent`` -- an UndirectedGraphModule

        - ``vector`` -- a dictionary, representing a sparse vector of coefficients with respect to the basis of ``parent``
        """
        self._parent = parent
        self._vector = defaultdict(lambda: self._parent.base_ring().zero())
        for key in vector:
            self._vector[key] = vector[key]

    def __repr__(self):
        """
        Return a string representation of ``self``.
        """
        terms = []
        for key in self._vector:
            c = self._vector[key]
            if c.is_zero():
                continue
            g, sign = self._parent._graph_basis.key_to_graph(key)
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
        Return the parent UndirectedGraphModule that ``self`` belongs to.
        """
        return self._parent

    def copy(self):
        """
        Return a copy of ``self``.
        """
        return self.__class__(self._parent, self._vector)

    __pos__ = copy

    def __neg__(self):
        """
        Return the negative of ``self``.
        """
        return self.__class__(self._parent, {k : -v for (k,v) in self._vector.items()})

    def __add__(self, other):
        """
        Return ``self`` added to ``other``.
        """
        if isinstance(other, self.__class__):
            v = self._vector.copy()
            for k in other._vector:
                v[k] += other._vector[k]
            return self.__class__(self._parent, v)
        elif other == 0:
            return self.copy()

    def __radd__(self, other):
        """
        Return ``other`` added to ``self``.
        """
        return self + other

    def __sub__(self, other):
        """
        Return ``other`` subtracted from ``self``.
        """
        if isinstance(other, self.__class__):
            v = self._vector.copy()
            for k in other._vector:
                v[k] -= other._vector[k]
            return self.__class__(self._parent, v)
        elif other == 0:
            return self.copy()

    def __rsub__(self, other):
        """
        Return ``self`` subtracted from ``other``.
        """
        return -(self - other)

    def __mul__(self, other):
        """
        Return ``self`` multiplied by ``other``.
        """
        v = self._vector.copy()
        if other in self._parent.base_ring():
            for k in v:
                v[k] *= other
            return self.__class__(self._parent, v)
        else:
            raise NotImplementedError

    def __rmul__(self, other):
        """
        Return ``other`` multiplied by ``self``.
        """
        return self * other

    def __eq__(self, other):
        """
        Return ``True`` if ``self`` is equal to ``other`` and ``False`` otherwise.
        """
        difference = self - other
        for k in difference._vector:
            if not difference._vector[k].is_zero():
                return False
        return True

    def bi_gradings(self):
        """
        Return the set of tuples ``(v,e)`` such that ``self`` contains terms with ``v`` vertices and ``e`` edges.
        """
        return set((key[0],key[1]) for key in self._vector)

    def homogeneous_part(self, vertices, edges):
        """
        Return the homogeneous part of ``self`` consisting only of terms with the given number of ``vertices`` and ``edges``.
        """
        v = {}
        for key in self._vector:
            if key[:2] == (vertices, edges):
                v[key] = self._vector[key]
        return self.__class__(self._parent, v)
    
    def insertion(self, position, other):
        """
        Return the insertion of ``other`` into ``self`` at the vertex ``position``.
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
                    term = UndirectedGraph(len(user) + len(victim) - 1, [tuple(e) for e in user_edges] + victim_edges)
                    terms.append([user_coeff*victim_coeff, term])
        return self._parent(terms)

class UndirectedGraphModule:
    """
    Module spanned by undirected graphs.
    """
    def __init__(self, base_ring, graph_basis):
        """
        Initialize ``self``.

        INPUT:

        - ``base_ring`` -- a ring, to be used as the ring of coefficients

        - ``graph_basis`` -- an UndirectedGraphBasis
        """
        self._base_ring = base_ring
        if not isinstance(graph_basis, UndirectedGraphBasis):
            raise ValueError('graph_basis must be UndirectedGraphBasis')
        self._graph_basis = graph_basis
        self.element_class = UndirectedGraphVector

    def base_ring(self):
        """
        Return the base ring of ``self``.
        """
        return self._base_ring

    def basis(self):
        """
        Return the basis of ``self``.
        """
        return self._graph_basis

    def __repr__(self):
        """
        Return a string representation of ``self``.
        """
        return 'Module over {} with {}'.format(self._base_ring, self._graph_basis)

    def __call__(self, arg):
        """
        Convert ``arg`` into an element of ``self``.
        """
        if isinstance(arg, UndirectedGraph):
            key, sign = self._graph_basis.graph_to_key(arg)
            if key is not None:
                return self.element_class(self, { key : self.base_ring().one() * sign })
            else: # must be zero
                return self.element_class(self, {})
        elif isinstance(arg, list):
            v = self.element_class(self, {})
            for term in arg:
                coeff, graph = term
                key, sign = self._graph_basis.graph_to_key(graph)
                coeff *= sign
                if key is not None:
                    v._vector[key] += coeff
            return v
        elif arg == 0:
            return self.element_class(self, {})
