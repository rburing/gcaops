from collections import defaultdict
from collections.abc import Iterable, MutableMapping
from itertools import combinations

class keydefaultdict(defaultdict):
    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        else:
            ret = self[key] = self.default_factory(key)
            return ret

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
        self._monomial_coefficients = monomial_coefficients

    def __repr__(self):
        """
        Return a string representation of ``self``.
        """
        terms = []
        for degree, coefficients in self._monomial_coefficients.items():
            for k, coefficient in enumerate(coefficients):
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

class SuperfunctionAlgebra:
    """
    Supercommutative algebra of superfunctions on a coordinate chart of a Z_2-graded space.

    Consisting of polynomials in the odd (degree 1) coordinates, with coefficients in the base ring (of even degree 0 functions).
    It is a free module over the base ring with an ordered basis consisting of sorted monomials in the odd coordinates.
    The elements encode skew-symmetric multi-derivations of the base ring, or multi-vectors.
    """
    def __init__(self, base_ring, even_coordinates=None, names=None):
        """
        Initialize ``self``.

        INPUT:

        - ``base_ring`` -- a commutative ring, considered as a ring of (even, degree 0) functions

        - ``even_coordinates`` -- (default: ``None``) a list or tuple of elements of ``base_ring``; if none is provided, then it is set to ``base_ring.gens()``

        - ``names`` -- (default: ``None``) a list or tuple of strings or a comma separated string, consisting of names for the odd coordinates; or a single string consisting of a prefix; if none is provided, then it is set to the prefix 'xi'
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
        self._gens = tuple(self.element_class(self, {1 : [1 if j == k else 0 for j in range(self.__ngens)]}) for k in range(self.__ngens))
        self._basis = keydefaultdict(lambda degree: list(combinations(range(self.__ngens), degree)))

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

    def _repr_monomial(self, degree, index):
        """
        Return a string representation of the respective monomial in the odd coordinates.

        INPUT:

        - ``degree`` -- a natural number, the degree of the monomial

        - ``index`` -- a natural number, the index of the monomial in the basis
        """
        return '*'.join(self._names[i] for i in self._basis[degree][index])
