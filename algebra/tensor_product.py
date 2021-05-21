from collections.abc import MutableSequence
from itertools import permutations
from util.permutation import selection_sort_graded
from math import factorial

class TensorProductElement:
    """
    Element of a tensor product of vector spaces.
    """
    def __init__(self, parent, terms):
        """
        Initialize this tensor product element.

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
        Return a string representation of this tensor product element.
        """
        return ' + '.join('⊗'.join('({})'.format(factor) for factor in term) for term in self._terms)

    def parent(self):
        """
        Return the parent TensorProduct that this tensor product element belongs to.
        """
        return self._parent

    def terms(self):
        """
        Return the list of terms of this tensor product element.
        """
        return self._terms

    def graded_symmetrization(self):
        """
        Return the graded symmetrization of this tensor product element.

        ASSUMPTION:

        Assumes each factor in each term of this tensor product element has a ``degree`` method and is homogeneous of that degree.
        """
        n = self._parent.nfactors()
        prefactor_inverse = factorial(n)
        new_terms = []
        for term in self._terms:
            symmetrized_terms = []
            # label the indices of identical factors by the minimal index which they are identical to
            minimal_label = {}
            for i in range(n):
                if not i in minimal_label:
                    for j in range(i,n):
                        if term[j] is term[i]:
                            minimal_label[j] = i
            # for each permutation, will use the minimal labeling as a key, to accumulate coefficients
            keys = []
            key_counter = 0
            accumulated_coeff = []
            degrees = [f.degree() for f in term]
            # sum over all permutations, but compute only once for each permutation of identical factors
            for sigma in permutations(range(n)):
                key = tuple(minimal_label[i] for i in sigma)
                sign = selection_sort_graded(list(sigma), degrees.copy())
                if not key in keys:
                    idx = key_counter
                    keys.append(key)
                    accumulated_coeff.append(0)
                    new_term = [term[k] for k in key]
                    new_term[0] /= prefactor_inverse
                    new_terms.append(new_term)
                    key_counter += 1
                else:
                    idx = keys.index(key)
                accumulated_coeff[idx] += sign
            # now distribute accumulated coefficients
            for idx in range(key_counter):
                new_terms[-key_counter+idx][0] *= accumulated_coeff[idx]
        return self.__class__(self._parent, new_terms)

class TensorProduct:
    """
    Tensor product of vector spaces.
    """
    def __init__(self, factors):
        """
        Initialize this tensor product.

        INPUT:

        - ``factors`` -- a list of vector spaces
        """
        self.element_class = TensorProductElement
        if not isinstance(factors, MutableSequence):
            raise ValueError('factors must be a list')
        self._factors = factors

    def __repr__(self):
        """
        Return a string representation of this tensor product.
        """
        return ' ⊗ '.join(repr(factor) for factor in self._factors)

    def factors(self):
        """
        Return the list of factors of this tensor product.
        """
        return self._factors

    def nfactors(self):
        """
        Return the number of factors of this tensor product.
        """
        return len(self._factors)

    def factor(self, index):
        """
        Return the ``index``th factor of this tensor product.
        """
        return self._factors[index]

    def __call__(self, arg):
        """
        Return ``arg`` converted into an element of this tensor product.
        """
        if isinstance(arg, self.element_class) and arg.parent() is self:
            return arg
        elif isinstance(arg, MutableSequence) and all(isinstance(term, MutableSequence) for term in arg):
            return self.element_class(self, arg)
        else:
            raise ValueError('cannot convert {} into element of {}'.format(arg, self))
