from util.permutation import selection_sort_graded
from abc import ABC, abstractmethod
from .tensor_product import TensorProduct
from functools import reduce
from itertools import combinations
import operator

# TODO: sum of operations

class SuperfunctionAlgebraOperation(ABC):
    """
    A homogeneous n-ary multi-linear operation acting on a SuperfunctionAlgebra.
    """
    def __init__(self, domain, codomain):
        """
        Initialize this operation.
        """
        self._domain = domain
        self._codomain = codomain

    def __repr__(self):
        """
        Return a string representation of this operation.
        """
        return 'Operation of arity {} and degree {} on {}'.format(self._domain.nfactors(), self.degree(), self._codomain)

    def domain(self):
        """
        Return the domain of this operation.
        """
        return self._domain

    def codomain(self):
        """
        Return the codomain of this operation.
        """
        return self._codomain

    @abstractmethod
    def degree(self):
        """
        Return the degree of this operation.
        """
        pass

    @abstractmethod
    def __call__(self, *arg):
        """
        Return the evaluation of this operation at ``arg``.
        """
        pass


class SuperfunctionAlgebraSymmetricOperation(SuperfunctionAlgebraOperation):
    """
    A homogeneous symmetric n-ary multi-linear operation acting on a SuperfunctionAlgebra.
    """
    def __repr__(self):
        """
        Return a string representation of this operation.
        """
        return 'Symmetric operation of arity {} and degree {} on {}'.format(self._domain.nfactors(), self.degree(), self._codomain)

    def bracket(self, other):
        """
        Return the Nijenhuis-Richardson bracket of this operation with the other operation.
        """
        return SuperfunctionAlgebraSymmetricBracketOperation(self, other)


class SuperfunctionAlgebraSymmetricBracketOperation(SuperfunctionAlgebraSymmetricOperation):
    """
    A homogeneous symmetric n-ary multi-linear operation acting on a SuperfunctionAlgebra, given by the Nijenhuis-Richardson bracket of two graded symmetric operations.
    """
    def __init__(self, *args):
        """
        Initialize this Nijenhuis-Richardson bracket.
        """
        if len(args) != 2:
            raise ValueError('The Nijenhuis-Richardson bracket takes two arguments.')
        if not all(isinstance(arg, SuperfunctionAlgebraSymmetricOperation) for arg in args):
            raise ValueError('The Nijenhuis-Richardson bracket is only defined for symmetric operations.')
        self.args = args
        self._codomain = args[0].codomain()
        superfunction_algebra = args[0].domain().factor(0)
        self._domain = superfunction_algebra.tensor_power(args[0].domain().nfactors() + args[1].domain().nfactors() - 1)

    def __repr__(self):
        """
        Return a string representation of this operation.
        """
        return 'Symmetric operation of arity {} and degree {} on {} given by the Nijenhuis-Richardson bracket of two symmetric operations'.format(self._domain.nfactors(), self.degree(), self._codomain)

    def degree(self):
        """
        Return the degree of this operation.
        """
        return self.args[0].degree() + self.args[1].degree()

    def __call__(self, *args):
        """
        Return the evaluation of this operation at the given arguments.
        """
        if len(args) == 1 and isinstance(args[0], self._domain.element_class) and args[0].parent() is self._domain:
            p = self.args[0].domain().nfactors() - 1
            q = self.args[1].domain().nfactors() - 1
            subtrahend_factor = 1 if (self.args[0].degree() % 2 == 1 and self.args[1].degree() % 2 == 1) else -1 # NOTE: shifted degrees
            result = self._codomain.zero()
            for term in args[0].terms(): # multi-linearity
                for sigma in combinations(range(p+q+1), q+1): # (q+1,p)-shuffles
                    sigma = sigma + tuple(k for k in range(p+q+1) if not k in sigma)
                    sign = selection_sort_graded(list(sigma), [v.degree() for v in term])
                    result += sign * self.args[0](*([self.args[1](*term[:q+1])] + term[q+1:]))
                for sigma in combinations(range(p+q+1), p+1): # (p+1,q)-shuffles
                    sigma = sigma + tuple(k for k in range(p+q+1) if not k in sigma)
                    sign = selection_sort_graded(list(sigma), [v.degree() for v in term])
                    result += sign * subtrahend_factor * self.args[1](*([self.args[0](*term[:p+1])] + term[p+1:]))
            return result
        elif len(args) == self._domain.nfactors():
            return self(self._domain([list(args)]))
        else:
            raise ValueError("input not recognized")


class SuperfunctionAlgebraSchoutenBracket(SuperfunctionAlgebraSymmetricOperation):
    """
    Schouten bracket on a SuperfunctionAlgebra.
    """
    def __init__(self, domain, codomain):
        """
        Initialize this Schouten bracket.
        """
        super().__init__(domain, codomain)

    def __repr__(self):
        """
        Return a string representation of this Schouten bracket.
        """
        return 'Schouten bracket on {}'.format(self._codomain)

    def degree(self):
        """
        Return the degree of this operation.
        """
        return -1

    def __call__(self, *arg):
        """
        Return the evaluation of this Schouten bracket at ``arg``.
        """
        if len(arg) == 2 and all(isinstance(arg[k], self._domain.factor(k).element_class) and arg[k].parent() is self._domain.factor(k) for k in range(2)):
            d = arg[0].degree()
            result = self._codomain.zero()
            for p in range(d+1):
                sign = -1 if p % 2 == 0 else 1 # NOTE: shifted degree
                result += sign * arg[0].homogeneous_part(p).bracket(arg[1])
            return result
        elif len(arg) == 1 and isinstance(arg[0], self._domain.element_class) and arg[0].parent() is self._domain:
            return sum(term[0].bracket(term[1]) for term in arg[0].terms())
        else:
            raise ValueError("input not recognized")


class SuperfunctionAlgebraUndirectedGraphOperation(SuperfunctionAlgebraOperation):
    """
    A homogeneous n-ary multi-linear operation acting on a SuperfunctionAlgebra, defined by a UndirectedGraphVector.
    """
    def __init__(self, domain, codomain, graph_vector):
        """
        Initialize this operation.
        """
        super().__init__(domain, codomain)
        arity = domain.nfactors()
        bi_gradings = graph_vector.bi_gradings()
        if len(bi_gradings) != 1:
            raise ValueError('graph_vector must be homogenous')
        bi_grading = next(iter(bi_gradings))
        if bi_grading[0] != arity:
            raise ValueError('graph_vector must have as many vertices as the number of factors in the domain')
        self._graph_vector = graph_vector
        self._degree = -bi_grading[1]

    def degree(self):
        """
        Return the degree of this operation.
        """
        return self._degree

    def _act_with_graph(self, graph, arg):
        """
        Return the evaluation of ``graph`` at ``arg``.

        ASSUMPTION:

        Assumes that each factor in each term of ``arg`` is homogeneous.
        """
        evens = self._codomain.even_coordinates()
        odds = self._codomain.odd_coordinates()
        terms = arg.terms()
        for e in graph.edges():
            new_terms = []
            for k in range(len(terms)):
                term0 = terms[k]
                if 0 in term0:
                    continue
                for k in range(self._codomain.ngens()):
                    left_odd_derivative = term0[e[1]].diff(odds[k])
                    if left_odd_derivative != 0:
                        left_even_derivative = term0[e[0]].diff(evens[k])
                        if left_even_derivative != 0:
                            left_term = [f.copy() for f in term0]
                            left_sign = 1 if sum(term0[j].degree() for j in range(e[1])) % 2 == 0 else -1
                            left_term[e[1]] = left_sign * left_odd_derivative
                            left_term[e[0]] = left_even_derivative
                            new_terms.append(left_term)
                    right_odd_derivative = term0[e[0]].diff(odds[k])
                    if right_odd_derivative != 0:
                        right_even_derivative = term0[e[1]].diff(evens[k])
                        if right_even_derivative != 0:
                            right_term = [f.copy() for f in term0]
                            right_sign = 1 if sum(term0[j].degree() for j in range(e[0])) % 2 == 0 else -1
                            right_term[e[0]] = right_sign * right_odd_derivative
                            right_term[e[1]] = right_even_derivative
                            new_terms.append(right_term)
            terms = new_terms
        return sum((reduce(operator.mul, term) for term in terms), self._codomain.zero())

    def __call__(self, *args):
        """
        Return the evaluation of this operation at ``args``.

        ASSUMPTION:

        Assumes that each factor in each term of ``args`` is homogeneous.
        """
        if len(args) == 1 and isinstance(args[0], self._domain.element_class) and args[0].parent() is self._domain:
            result = self._codomain.zero()
            for (c, g) in self._graph_vector:
                result += c*self._act_with_graph(g, args[0])
            return result
        elif len(args) == self._domain.nfactors():
            return self(self._domain([list(args)]))
        else:
            raise ValueError("input not recognized")


class SuperfunctionAlgebraSymmetricUndirectedGraphOperation(SuperfunctionAlgebraUndirectedGraphOperation, SuperfunctionAlgebraSymmetricOperation):
    """
    A homogeneous n-ary multi-linear symmetric operation acting on a SuperfunctionAlgebra, defined by a UndirectedGraphVector.
    """
    def __call__(self, *args):
        """
        Return the evaluation of this operation at `args``.

        ASSUMPTION:

        Assumes that each factor in each term of ``args`` is homogeneous.
        """
        if len(args) == 1 and isinstance(args[0], self._domain.element_class) and args[0].parent() is self._domain:
            return super().__call__(args[0].graded_symmetrization())
        elif len(args) == self._domain.nfactors():
            return self(self._domain([list(args)]))
        else:
            raise ValueError("input not recognized")


class SuperfunctionAlgebraDirectedGraphOperation(SuperfunctionAlgebraOperation):
    """
    A homogeneous n-ary multi-linear operation on a SuperfunctionAlgebra, defined by a DirectedGraphVector.
    """
    def __init__(self, domain, codomain, graph_vector):
        """
        Initialize this operation.
        """
        super().__init__(domain, codomain)
        arity = domain.nfactors()
        bi_gradings = graph_vector.bi_gradings()
        if len(bi_gradings) != 1:
            raise ValueError('graph_vector must be homogenous')
        bi_grading = next(iter(bi_gradings))
        if bi_grading[0] != arity:
            raise ValueError('graph_vector must have as many vertices as the number of factors in the domain')
        self._graph_vector = graph_vector
        self._degree = -bi_grading[1]

    def degree(self):
        """
        Return the degree of this operation.
        """
        return self._degree

    def _act_with_graph(self, graph, arg):
        """
        Return the evaluation of ``graph`` at ``arg``.

        ASSUMPTION:

        Assumes that each factor in each term of ``arg`` is homogeneous.
        """
        evens = self._codomain.even_coordinates()
        odds = self._codomain.odd_coordinates()
        terms = arg.terms()
        for e in graph.edges():
            new_terms = []
            for k in range(len(terms)):
                term0 = terms[k]
                if 0 in term0:
                    continue
                for k in range(self._codomain.ngens()):
                    right_odd_derivative = term0[e[0]].diff(odds[k])
                    if right_odd_derivative != 0:
                        right_even_derivative = term0[e[1]].diff(evens[k])
                        if right_even_derivative != 0:
                            right_term = [f.copy() for f in term0]
                            right_sign = 1 if sum(term0[j].degree() for j in range(e[0])) % 2 == 0 else -1
                            right_term[e[0]] = right_sign * right_odd_derivative
                            right_term[e[1]] = right_even_derivative
                            new_terms.append(right_term)
            terms = new_terms
        return sum((reduce(operator.mul, term) for term in terms), self._codomain.zero())

    def __call__(self, *args):
        """
        Return the evaluation of this operation at ``args``.

        ASSUMPTION:

        Assumes that each factor in each term of ``args`` is homogeneous.
        """
        if len(args) == 1 and isinstance(args[0], self._domain.element_class) and args[0].parent() is self._domain:
            result = self._codomain.zero()
            for (c, g) in self._graph_vector:
                result += c*self._act_with_graph(g, args[0])
            return result
        elif len(args) == self._domain.nfactors():
            return self(self._domain([list(args)]))
        else:
            raise ValueError("input not recognized")


class SuperfunctionAlgebraSymmetricDirectedGraphOperation(SuperfunctionAlgebraDirectedGraphOperation, SuperfunctionAlgebraSymmetricOperation):
    """
    A homogeneous symmetric n-ary multi-linear operation acting on a SuperfunctionAlgebra, defined by a DirectedGraphVector.
    """
    def __call__(self, *args):
        """
        Return the evaluation of this operation at `args``.

        ASSUMPTION:

        Assumes that each factor in each term of ``args`` is homogeneous.
        """
        if len(args) == 1 and isinstance(args[0], self._domain.element_class) and args[0].parent() is self._domain:
            return super().__call__(args[0].graded_symmetrization())
        elif len(args) == self._domain.nfactors():
            return self(self._domain([list(args)]))
        else:
            raise ValueError("input not recognized")
