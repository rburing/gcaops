from collections.abc import Iterable, MutableMapping
from itertools import combinations
from util.permutation import selection_sort
from util.misc import keydefaultdict
from .superfunction_algebra_operation import SuperfunctionAlgebraSchoutenBracket
from .tensor_product import TensorProduct

class Superfunction:
    """
    Superfunction on a coordinate chart of a Z_2-graded space.

    A polynomial in the odd coordinates, with coefficients in the base ring (of even degree 0 functions).
    """
    def __init__(self, parent, monomial_coefficients):
        """
        Initialize ``self``.

        INPUT:

        - ``parent`` - a SuperfunctionAlgebra (which has an ordered basis of monomials in the odd coordinates)

        - ``monomial_coefficients`` - a dictionary, taking a natural number ``d`` to a list of coefficients of the monomials of degree ``d`` in the ordered basis of ``parent``
        """
        if not isinstance(parent, SuperfunctionAlgebra):
            raise TypeError('parent must be a SuperfunctionAlgebra')
        self._parent = parent
        if not isinstance(monomial_coefficients, MutableMapping):
            raise TypeError('monomial_coefficients must be a dictionary')
        self._monomial_coefficients = keydefaultdict(lambda degree: [self._parent.base_ring().zero() for k in range(self._parent.dimension(degree))])
        for degree in monomial_coefficients:
            self._monomial_coefficients[degree] = monomial_coefficients[degree]
            for k in range(len(self._monomial_coefficients[degree])):
                self._monomial_coefficients[degree][k] = self._parent.base_ring()(self._monomial_coefficients[degree][k]) # conversion

    def __repr__(self):
        """
        Return a string representation of ``self``.
        """
        terms = []
        for degree in reversed(sorted(self._monomial_coefficients.keys())):
            for k, coefficient in enumerate(self._monomial_coefficients[degree]):
                c = repr(coefficient)
                if c == '0':
                    continue
                elif c == '1' and degree > 0: # mainly for generators and basis
                    term = self._parent._repr_monomial(degree, k)
                elif degree == 0:
                    term = '({})'.format(c)
                else:
                    term = '({})*{}'.format(c, self._parent._repr_monomial(degree, k))
                terms.append(term)
        if terms:
            return ' + '.join(terms)
        else:
            return '0'

    def parent(self):
        """
        Return the parent SuperfunctionAlgebra that ``self`` belongs to.
        """
        return self._parent

    def homogeneous_part(self, degree):
        """
        Return the homogeneous part of ``self`` of total degree ``degree`` in the odd coordinates.

        NOTE::

            Returns a Superfunction whose homogeneous component of degree ``degree`` is a *reference* to the respective component of ``self``.
        """
        return self.__class__(self._parent, { degree : self._monomial_coefficients[degree] })

    def map_coefficients(self, f):
        """
        """
        monomial_coefficients = keydefaultdict(lambda degree: [self._parent.base_ring().zero() for k in range(self._parent.dimension(degree))])
        for degree in self._monomial_coefficients:
            for k in range(len(self._monomial_coefficients[degree])):
                monomial_coefficients[degree][k] = self._parent._simplify(f(self._monomial_coefficients[degree][k]))
        return self.__class__(self._parent, monomial_coefficients)

    def copy(self):
        """
        Return a copy of ``self``.
        """
        return self.map_coefficients(lambda c: c)

    __pos__ = copy

    def __neg__(self):
        """
        Return the negative of ``self``.
        """
        return self.map_coefficients(lambda c: -c)

    def __add__(self, other):
        """
        Return ``self`` added to ``other``.
        """
        monomial_coefficients = keydefaultdict(lambda degree: [self._parent.base_ring().zero() for k in range(self._parent.dimension(degree))])
        for degree in self._monomial_coefficients:
            for k in range(len(self._monomial_coefficients[degree])):
                monomial_coefficients[degree][k] = self._monomial_coefficients[degree][k]
        if isinstance(other, self.__class__):
            for degree in other._monomial_coefficients:
                for k in range(len(other._monomial_coefficients[degree])):
                    monomial_coefficients[degree][k] = self._parent._simplify(monomial_coefficients[degree][k] + other._monomial_coefficients[degree][k])
        elif other in self._parent.base_ring():
            monomial_coefficients[0][0] = self._parent._simplify(monomial_coefficients[0][0] + other)
        else:
            raise NotImplementedError
        return self.__class__(self._parent, monomial_coefficients)

    def __radd__(self, other):
        """
        Return ``other`` added to ``self``.
        """
        return self + other

    def __sub__(self, other):
        """
        Return ``self`` minus ``other``.
        """
        monomial_coefficients = keydefaultdict(lambda degree: [self._parent.base_ring().zero() for k in range(self._parent.dimension(degree))])
        for degree in self._monomial_coefficients:
            for k in range(len(self._monomial_coefficients[degree])):
                monomial_coefficients[degree][k] = self._monomial_coefficients[degree][k]
        if isinstance(other, self.__class__):
            for degree in other._monomial_coefficients:
                for k in range(len(other._monomial_coefficients[degree])):
                    monomial_coefficients[degree][k] = self._parent._simplify(monomial_coefficients[degree][k] - other._monomial_coefficients[degree][k])
        elif other in self._parent.base_ring():
            monomial_coefficients[0][0] -= other
        else:
            raise NotImplementedError
        return self.__class__(self._parent, monomial_coefficients)

    def __rsub__(self, other):
        """
        Return ``other`` minus ``self``.
        """
        return -(self - other)

    def __mul__(self, other):
        """
        Return ``self`` multiplied by ``other``.
        """
        monomial_coefficients = keydefaultdict(lambda degree: [self._parent.base_ring().zero() for k in range(self._parent.dimension(degree))])
        if isinstance(other, self.__class__):
            for degree1 in self._monomial_coefficients:
                for k1 in range(len(self._monomial_coefficients[degree1])):
                    if self._parent._is_zero(self._monomial_coefficients[degree1][k1]):
                        continue
                    for degree2 in other._monomial_coefficients:
                        for k2 in range(len(other._monomial_coefficients[degree2])):
                            if self._parent._is_zero(other._monomial_coefficients[degree2][k2]):
                                continue
                            prod, sign = self._parent._mul_on_basis(degree1,k1,degree2,k2)
                            if prod is not None:
                                monomial_coefficients[degree1+degree2][prod] = self._parent._simplify(monomial_coefficients[degree1+degree2][prod] + sign * self._monomial_coefficients[degree1][k1] * other._monomial_coefficients[degree2][k2])
        elif other in self._parent.base_ring():
            for degree in self._monomial_coefficients:
                for k in range(len(self._monomial_coefficients[degree])):
                    monomial_coefficients[degree][k] = self._parent._simplify(self._monomial_coefficients[degree][k] * other)
        else:
            raise NotImplementedError
        return self.__class__(self._parent, monomial_coefficients)

    def __rmul__(self, other):
        """
        Return ``other`` multiplied by ``self``.

        NOTE::

            This assumes that ``other`` commutes with ``self``.
            It is justified because this function only gets called when ``other`` is even.
        """
        return self * other

    def __truediv__(self, other):
        """
        Return ``self`` divided by ``other``.
        """
        return self.map_coefficients(lambda c: c / other)

    def __eq__(self, other):
        """
        Return ``True`` if ``self`` equals ``other`` and ``False`` otherwise.
        """
        difference = self - other
        for degree in difference._monomial_coefficients:
            for k in range(len(difference._monomial_coefficients[degree])):
                if not self._parent._is_zero(difference._monomial_coefficients[degree][k]):
                    return False
        return True

    def degree(self):
        """
        Return the degree of ``self`` as a polynomial in the odd coordinates.
        """
        for d in reversed(sorted(self._monomial_coefficients.keys())):
            for k in range(len(self._monomial_coefficients[d])):
                if any(not self._parent._is_zero(c) for c in self._monomial_coefficients[d]):
                    return d
        return 0

    def derivative(self, *args):
        """
        Return the derivative of ``self`` with respect to ``args``.

        INPUT:

        - ``args`` -- an odd coordinate or an even coordinate, or a list of such
        """
        if len(args) > 1:
            result = self
            for arg in args:
                result = result.derivative(arg)
            return result
        elif len(args) == 1 and any(args[0] is xi for xi in self._parent.gens()):
            j = self._parent.gens().index(args[0])
            monomial_coefficients = keydefaultdict(lambda degree: [self._parent.base_ring().zero() for k in range(self._parent.dimension(degree))])
            for degree in self._monomial_coefficients:
                for k in range(len(self._monomial_coefficients[degree])):
                    derivative, sign = self._parent._derivative_on_basis(degree, k, j)
                    if derivative is not None:
                        monomial_coefficients[degree-1][derivative] = self._parent._simplify(sign * self._monomial_coefficients[degree][k])
            return self.__class__(self._parent, monomial_coefficients)
        elif len(args) == 1 and any(args[0] is x for x in self._parent.even_coordinates()):
            monomial_coefficients = keydefaultdict(lambda degree: [self._parent.base_ring().zero() for k in range(self._parent.dimension(degree))])
            for degree in self._monomial_coefficients:
                for k in range(len(self._monomial_coefficients[degree])):
                    monomial_coefficients[degree][k] = self._parent._simplify(self._monomial_coefficients[degree][k].derivative(args[0]))
            return self.__class__(self._parent, monomial_coefficients)
        elif len(args) == 1:
            # by now we know args[0] is not identically a coordinate, but maybe it is equal to one:
            try:
                actual_xi_idx = self._parent.gens().index(args[0])
                return self.derivative(self._parent.gen(actual_xi_idx))
            except ValueError:
                try:
                    actual_x_idx = self._parent.even_coordinates().index(args[0])
                    return self.derivative(self._parent.even_coordinate(actual_x_idx))
                except ValueError:
                    raise ValueError("{} not recognized as a coordinate".format(args[0]))
        else:
            raise ValueError("Don't know how to take derivative with respect to {}".format(args))

    diff = derivative

    def bracket(self, other):
        """
        Return the Schouten bracket (odd Poisson bracket) of ``self`` with ``other``.
        """
        if len(self._monomial_coefficients.keys()) == 1: # first argument is homogeneous
            degree = next(iter(self._monomial_coefficients))
            sign = 1 if degree % 2 == 0 else -1
            return sum(self.derivative(self._parent.even_coordinate(i))*other.derivative(self._parent.gen(i)) + sign*self.derivative(self._parent.gen(i))*other.derivative(self._parent.even_coordinate(i)) for i in range(self._parent.ngens()))
        else:
            return sum(self.homogeneous_part(d).bracket(other) for d in self._monomial_coefficients.keys())

class SuperfunctionAlgebra:
    """
    Supercommutative algebra of superfunctions on a coordinate chart of a Z_2-graded space.

    Consisting of polynomials in the odd (degree 1) coordinates, with coefficients in the base ring (of even degree 0 functions).
    It is a free module over the base ring with an ordered basis consisting of sorted monomials in the odd coordinates.
    The elements encode skew-symmetric multi-derivations of the base ring, or multi-vectors.
    """
    def __init__(self, base_ring, even_coordinates=None, names=None, simplify=lambda x: x, is_zero=lambda x: x.is_zero()):
        """
        Initialize ``self``.

        INPUT:

        - ``base_ring`` -- a commutative ring, considered as a ring of (even, degree 0) functions

        - ``even_coordinates`` -- (default: ``None``) a list or tuple of elements of ``base_ring``; if none is provided, then it is set to ``base_ring.gens()``

        - ``names`` -- (default: ``None``) a list or tuple of strings or a comma separated string, consisting of names for the odd coordinates; or a single string consisting of a prefix; if none is provided, then it is set to the prefix 'xi'

        - ``simplify`` -- (default: ``lambda x: x``) a function that returns a simplification of an element of the base ring (will be used in each operation on elements that affects coefficients)

        - ``is_zero`` -- (default: ``lambda x: x.is_zero()``) a function that returns ``True`` when a simplified element of the base ring is equal to zero (will be used to decide equality of elements, to calculate the degree of elements, and to skip terms in some operations on elements)
        """
        self.element_class = Superfunction
        self._base_ring = base_ring
        if even_coordinates:
            self._even_coordinates = even_coordinates
        elif hasattr(base_ring, 'gens'):
            self._even_coordinates = base_ring.gens()
        else:
            raise ValueError('Even coordinates not specified and could not be determined from base ring')
        if names is None:
            names = 'xi'
        if isinstance(names, str):
            if ',' in names:
                names = names.split(',')
            else:
                names = ['{}{}'.format(names, k) for k in range(len(self._even_coordinates))]
        elif not isinstance(names, Iterable) or not all(isinstance(name, str) for name in names):
            raise ValueError('Format of odd coordinate names {} not recognized'.format(names))
        if len(names) != len(self._even_coordinates):
            raise ValueError("Number of odd coordinate names in {} does not match number of even coordinates".format(names))
        self._names = tuple(names)
        self.__ngens = len(names)
        self._gens = tuple(self.element_class(self, {1 : [self._base_ring.one() if j == k else self._base_ring.zero() for j in range(self.__ngens)]}) for k in range(self.__ngens))
        self._basis = keydefaultdict(lambda degree: list(combinations(range(self.__ngens), degree)))
        self._simplify = simplify
        self._is_zero = is_zero
        self._tensor_powers = keydefaultdict(lambda n: TensorProduct([self]*n))
        self._schouten_bracket = SuperfunctionAlgebraSchoutenBracket(self._tensor_powers[2], self)

    def __repr__(self):
        """
        Return a string representation of ``self``.
        """
        return "Superfunction algebra over {} with even coordinates {} and odd coordinates {}".format(self._base_ring, self._even_coordinates, self._gens)

    def __call__(self, arg):
        """
        Return ``arg`` converted into an element of ``self``.
        """
        if arg in self._base_ring:
            return self.element_class(self, { 0 : [arg]})
        elif isinstance(arg, self.element_class) and arg.parent() == self:
            return arg
        else:
            raise ValueError('cannot convert {} into element of {}'.format(arg, self))

    def base_ring(self):
        """
        Return the base ring of ``self``, consisting of (even, degree 0) functions.
        """
        return self._base_ring

    def even_coordinates(self):
        """
        Return the even coordinates in the base ring of ``self``.
        """
        return self._even_coordinates

    def even_coordinate(self, i):
        """
        Return the ``i``th even coordinate in the base ring of ``self``.
        """
        return self._even_coordinates[i]

    def ngens(self):
        """
        Return the number of odd coordinates of ``self``.
        """
        return self.__ngens

    def _first_ngens(self, n):
        """
        Return the first ``n`` odd coordinates of ``self``.
        """
        return self._gens[:n]

    def gens(self):
        """
        Return the tuple of odd coordinates of ``self``.
        """
        return self._gens

    odd_coordinates = gens

    def gen(self, i):
        """
        Return the ``i``th odd coordinate of ``self``.
        """
        return self._gens[i]

    odd_coordinate = gen

    def _repr_monomial(self, degree, index):
        """
        Return a string representation of the respective monomial in the odd coordinates.

        INPUT:

        - ``degree`` -- a natural number, the degree of the monomial

        - ``index`` -- a natural number, the index of the monomial in the basis
        """
        return '*'.join(self._names[i] for i in self._basis[degree][index])

    def dimension(self, degree):
        """
        Return the dimension of the graded component spanned by monomials of the given ``degree`` in the odd coordinates (as a module over the base ring).

        INPUT:

        - ``degree`` -- a natural number
        """
        return len(self._basis[degree])

    def _mul_on_basis(self, degree1, k1, degree2, k2):
        """
        Return the index and the sign of the monomial that results from multiplying the ``k1``th monomial of degree ``degree1`` by the ``k2``th monomial of degree ``degree2``.
        """
        if degree1 + degree2 > self.__ngens:
            return None, 1
        left = self._basis[degree1][k1]
        right = self._basis[degree2][k2]
        lst = list(left+right)
        sign = selection_sort(lst)
        # detect repetitions in sorted list
        for i in range(0, len(lst)-1):
            if lst[i] == lst[i+1]:
                return None, 1
        prod = tuple(lst)
        assert prod in self._basis[degree1+degree2]
        return self._basis[degree1+degree2].index(prod), sign

    def _derivative_on_basis(self, degree, i, j):
        """
        Return the index and the sign of the derivative of the ``i``th monomial of degree ``degree`` in the basis, with respect to the ``j``th odd coordinate.
        """
        monomial = self._basis[degree][i]
        if not j in monomial:
            return None, 1
        lst = list(monomial)
        sign = 1 if lst.index(j) % 2 == 0 else -1
        lst.remove(j) # remove first instance
        sign *= selection_sort(lst)
        derivative = tuple(lst)
        assert derivative in self._basis[degree-1]
        return self._basis[degree-1].index(derivative), sign

    def tensor_power(self, n):
        """
        Return the ``n``th tensor power of ``self``.
        """
        return self._tensor_powers[n]

    def schouten_bracket(self):
        """
        Return the Schouten bracket on ``self``.
        """
        return self._schouten_bracket

