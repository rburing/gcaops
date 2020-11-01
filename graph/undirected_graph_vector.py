from .undirected_graph import UndirectedGraph
from .undirected_graph_basis import UndirectedGraphBasis
from collections import defaultdict

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
        for idx in self._vector:
            c = self._vector[idx]
            if c.is_zero():
                continue
            g, sign = self._parent._graph_basis.index_to_graph(idx)
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
        return __class__(self._parent, self._vector)

    __pos__ = copy

    def __neg__(self):
        """
        Return the negative of ``self``.
        """
        return __class__(self._parent, {k : -v for (k,v) in self._vector.items()})

    def __add__(self, other):
        """
        Return ``self`` added to ``other``.
        """
        v = self._vector.copy()
        for k in other._vector:
            v[k] += other._vector[k]
        return __class__(self._parent, v)

    def __sub__(self, other):
        """
        Return ``other`` subtracted from ``self``.
        """
        v = self._vector.copy()
        for k in other._vector:
            v[k] -= other._vector[k]
        return __class__(self._parent, v)

    def __mul__(self, other):
        """
        Return ``self`` multiplied by ``other``.
        """
        v = self._vector.copy()
        if other in self._parent.base_ring():
            for k in v:
                v[k] *= other
            return __class__(self._parent, v)
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
            index, sign = self._graph_basis.graph_to_index(arg)
            if index is not None:
                return self.element_class(self, { index : self.base_ring().one() * sign })
            else: # must be zero
                return self.element_class(self, {})
