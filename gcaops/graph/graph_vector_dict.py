r"""
Graph vector (dictionary backend)
"""
from collections import defaultdict
from functools import partial
from .graph_vector import GraphVector, GraphModule
from .graph_basis import GraphBasis

def zero_in_base_ring(graph_module):
    return graph_module.base_ring().zero()

class GraphVector_dict(GraphVector):
    """
    Vector representing a linear combination of graphs (stored as a dictionary).
    """
    def __init__(self, parent, vector):
        """
        Initialize this graph vector.

        INPUT:

        - ``parent`` -- a :class:`~gcaops.graph.graph_vector.GraphModule`

        - ``vector`` -- a dictionary, representing a sparse vector of coefficients with respect to the basis of ``parent``
        """
        assert isinstance(parent, GraphModule_dict)
        self._parent = parent
        self._vector = defaultdict(partial(zero_in_base_ring, self._parent))
        for key in vector:
            self._vector[key] = vector[key]

    def __repr__(self):
        """
        Return a string representation of this graph vector.
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
        Return the parent :class:`~gcaops.graph.graph_vector.GraphModule` that this graph vector belongs to.
        """
        return self._parent

    def copy(self):
        """
        Return a copy of this graph vector.
        """
        return self.__class__(self._parent, self._vector)

    def __iter__(self):
        """
        Facilitates iterating over this graph vector, yielding tuples of the form ``(coeff, graph)``.
        """
        for key in self._vector:
            c = self._vector[key]
            if c.is_zero():
                continue
            g, sign = self._parent._graph_basis.key_to_graph(key)
            c *= sign
            yield (c, g)

    def __len__(self):
        """
        Return the number of graphs with nonzero coefficients in this graph vector.
        """
        count = 0
        for c in self._vector.values():
            if not c.is_zero():
                count += 1
        return count

    __pos__ = copy

    def __neg__(self):
        """
        Return the negative of this graph vector.
        """
        return self.__class__(self._parent, {k : -v for (k,v) in self._vector.items()})

    def __add__(self, other):
        """
        Return this graph vector added to ``other``.
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
        Return ``other`` added to this graph vector.
        """
        return self + other

    def __sub__(self, other):
        """
        Return ``other`` subtracted from this graph vector.
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
        Return this graph vector subtracted from ``other``.
        """
        return -(self - other)

    def __mul__(self, other):
        """
        Return this graph vector multiplied by ``other``.
        """
        if other in self._parent.base_ring():
            v = self._vector.copy()
            for k in v:
                v[k] *= other
            return self.__class__(self._parent, v)
        else:
            return NotImplemented

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
        for k in difference._vector:
            if not difference._vector[k].is_zero():
                return False
        return True

    def gradings(self):
        """
        Return the set of grading tuples such that this graph vector contains terms with those gradings.
        """
        grading_size = self._parent.basis().grading_size
        return set(key[:grading_size] for key in self._vector)

    def homogeneous_part(self, *grading):
        """
        Return the homogeneous part of this graph vector consisting only of terms with the given ``grading``.
        """
        grading_size = self._parent.basis().grading_size
        v = {}
        for key in self._vector:
            if key[:grading_size] == grading:
                v[key] = self._vector[key]
        return self.__class__(self._parent, v)

    def map_coefficients(self, f, new_parent=None):
        """
        Apply ``f`` to each of this graph vector's coefficients and return the resulting graph vector.
        """
        if new_parent is None:
            new_parent = self._parent
        new_vector = {}
        for k,v in self._vector.items():
            new_v = f(v)
            if new_v.is_zero():
                continue
            new_vector[k] = new_v
        return self.__class__(new_parent, new_vector)

    def insertion(self, position, other, **kwargs):
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
                product_coeff = user_coeff * victim_coeff
                if product_coeff.is_zero():
                    continue
                user, user_sign = self._parent._graph_basis.key_to_graph(user_key)
                user_coeff *= user_sign
                victim, victim_sign = other._parent._graph_basis.key_to_graph(victim_key)
                victim_coeff *= victim_sign
                for g in user._insertion_graphs(position, victim, **kwargs):
                    terms.append([product_coeff, g])
        return self._parent(terms)

class GraphModule_dict(GraphModule):
    """
    Module spanned by graphs (with elements stored as dictionaries).
    """
    def __init__(self, base_ring, graph_basis):
        """
        Initialize this graph module.

        INPUT:

        - ``base_ring`` -- a ring, to be used as the ring of coefficients

        - ``graph_basis`` -- a :class:`~gcaops.graph.graph_basis.GraphBasis`
        """
        self._base_ring = base_ring
        if not isinstance(graph_basis, GraphBasis):
            raise ValueError('graph_basis must be a GraphBasis')
        self._graph_basis = graph_basis
        self.element_class = GraphVector_dict

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
        if isinstance(arg, self._graph_basis.graph_class):
            key, sign = self._graph_basis.graph_to_key(arg)
            if key is not None:
                return self.element_class(self, { key : self.base_ring().one() * sign })
            else: # must be zero
                return self.zero()
        elif isinstance(arg, list):
            v = self.zero()
            for term in arg:
                coeff, graph = term
                key, sign = self._graph_basis.graph_to_key(graph)
                coeff *= sign
                if key is not None:
                    v._vector[key] += self._base_ring(coeff)
            return v
        elif isinstance(arg, self.element_class) and arg.parent() is self:
            return arg
        elif isinstance(arg, GraphVector):
            return self(list(arg))
        elif arg == 0:
            return self.zero()
        else:
            raise NotImplementedError
