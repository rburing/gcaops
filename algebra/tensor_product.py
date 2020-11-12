from collections.abc import MutableSequence
from itertools import permutations
from util.permutation import selection_sort_graded

class TensorProductElement:
    """
    Element of a tensor product of vector spaces.
    """
    def __init__(self, parent, terms):
        """
        Initialize ``self``.

        INPUT:

        - ``terms`` - list of lists
        """
        self._parent = parent
        if isinstance(terms, MutableSequence) and all(isinstance(term, MutableSequence) for term in terms):
            if all(len(term) == self._parent.nfactors() for term in terms):
                self._terms = terms
                for n in range(len(terms)):
                    for m in range(len(terms[n])):
                        self._terms[n][m] = self._parent.factor(m)(terms[n][m]) # conversion
            else:
                raise ValueError('each list in terms must be of length {}'.format(self._parent.nfactors()))
        else:
            raise ValueError('terms must be a list of lists')

    def __repr__(self):
        """
        Return a string representation of ``self``.
        """
        return ' + '.join('⊗'.join('({})'.format(factor) for factor in term) for term in self._terms)

    def parent(self):
        """
        Return the parent TensorProduct that ``self`` belongs to.
        """
        return self._parent

    def terms(self):
        """
        Return the list of terms of ``self``.
        """
        return self._terms

    def graded_symmetrization(self):
        """
        Return the graded symmetrization of ``self``.

        ASSUMPTION:

        Assumes each factor in each term of ``self`` has a ``degree`` method and is homogeneous of that degree.
        """
        # TODO: optimize
        n = self._parent.nfactors()
        new_terms = []
        for term in self._terms:
            degrees = [f.degree() for f in term]
            for sigma in permutations(range(n)):
                new_term = [term[sigma[k]] for k in range(n)]
                sign = selection_sort_graded(list(sigma), degrees.copy())
                new_term[0] *= sign
                new_terms.append(new_term)
        return self.__class__(self._parent, new_terms)

class TensorProduct:
    """
    Tensor product of vector spaces.
    """
    def __init__(self, factors):
        """
        Initialize ``self``.

        INPUT:

        - ``factors`` -- a list of vector spaces
        """
        self.element_class = TensorProductElement
        if not isinstance(factors, MutableSequence):
            raise ValueError('factors must be a list')
        self._factors = factors

    def __repr__(self):
        """
        Return a string representation of ``self``.
        """
        return ' ⊗ '.join(repr(factor) for factor in self._factors)

    def factors(self):
        """
        Return the list of factors of ``self``.
        """
        return self._factors

    def nfactors(self):
        """
        Return the number of factors of ``self``.
        """
        return len(self._factors)

    def factor(self, index):
        """
        Return the ``index``th factor of ``self``.
        """
        return self._factors[index]

    def __call__(self, arg):
        """
        Return ``arg`` converted into an element of ``self``.
        """
        if isinstance(arg, self.element_class) and arg.parent() == self:
            return arg
        elif isinstance(arg, MutableSequence) and all(isinstance(term, MutableSequence) for term in arg):
            return self.element_class(self, arg)
        else:
            raise ValueError('cannot convert {} into element of {}'.format(arg, self))
