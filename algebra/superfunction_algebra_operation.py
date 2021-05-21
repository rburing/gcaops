from abc import ABC, abstractmethod
from .tensor_product import TensorProduct
from functools import reduce
import operator

# TODO: bracket on the space of operations, sum of operations.

class SuperfunctionAlgebraOperation(ABC):
    """
    An n-ary multi-linear operation acting on a SuperfunctionAlgebra.
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
        return 'Operation of arity {} on {}'.format(self._domain.nfactors(), self._codomain)

    @abstractmethod
    def __call__(self, *arg):
        """
        Return the evaluation of this operation at ``arg``.
        """
        pass

class SuperfunctionAlgebraSchoutenBracket(SuperfunctionAlgebraOperation):
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

    def __call__(self, *arg):
        """
        Return the evaluation of this Schouten bracket at ``arg``.
        """
        if len(arg) == 2 and all(isinstance(arg[k], self._domain.factor(k).element_class) and arg[k].parent() is self._domain.factor(k) for k in range(2)):
            return arg[0].bracket(arg[1])
        elif len(arg) == 1 and isinstance(arg[0], self._domain.element_class) and arg[0].parent() is self._domain:
            return sum(term[0].bracket(term[1]) for term in arg[0].terms())
        else:
            raise ValueError("input not recognized")

class SuperfunctionAlgebraUndirectedGraphOperation(SuperfunctionAlgebraOperation):
    """
    Operation on a SuperfunctionAlgebra defined by a UndirectedGraphVector.
    """
    def __init__(self, domain, codomain, graph_vector):
        """
        Initialize this operation.
        """
        super().__init__(domain, codomain)
        self._graph_vector = graph_vector

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

class SuperfunctionAlgebraSymmetricUndirectedGraphOperation(SuperfunctionAlgebraUndirectedGraphOperation):
    """
    Graded symmetric operation on a SuperfunctionAlgebra defined by a UndirectedGraphVector.
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
    Operation on a SuperfunctionAlgebra defined by a DirectedGraphVector.
    """
    def __init__(self, domain, codomain, graph_vector):
        """
        Initialize this operation.
        """
        super().__init__(domain, codomain)
        self._graph_vector = graph_vector

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

class SuperfunctionAlgebraSymmetricDirectedGraphOperation(SuperfunctionAlgebraDirectedGraphOperation):
    """
    Graded symmetric operation on a SuperfunctionAlgebra defined by a DirectedGraphVector.
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
