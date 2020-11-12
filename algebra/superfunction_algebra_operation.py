from abc import ABC, abstractmethod
from .tensor_product import TensorProduct
from functools import reduce
import operator

class SuperfunctionAlgebraOperation(ABC):
    """
    An n-ary multi-linear operation acting on a SuperfunctionAlgebra.
    """
    def __init__(self, domain, codomain):
        """
        Initialize ``self``.
        """
        self._domain = domain
        self._codomain = codomain

    def __repr__(self):
        """
        Return a string representation of ``self``.
        """
        return 'Operation of arity {} on {}'.format(self._domain.nfactors(), self._codomain)

    @abstractmethod
    def __call__(self, *arg):
        """
        Return the evaluation of ``self`` at ``arg``.
        """
        pass

class SuperfunctionAlgebraSchoutenBracket(SuperfunctionAlgebraOperation):
    """
    Schouten bracket on a SuperfunctionAlgebra.
    """
    def __init__(self, domain, codomain):
        """
        Initialize ``self``.
        """
        super().__init__(domain, codomain)

    def __repr__(self):
        """
        Return a string representation of ``self``.
        """
        return 'Schouten bracket on {}'.format(self._codomain)

    def __call__(self, *arg):
        """
        Return the evaluation of ``self`` at ``arg``.
        """
        if len(arg) == 2 and all(isinstance(arg[k], self._domain.factor(k).element_class) and arg[k].parent() is self._domain.factor(k) for k in range(2)):
            return arg[0].bracket(arg[1])
        elif len(arg) == 1 and isinstance(arg[0], self._domain.element_class) and arg[0].parent() is self._domain:
            return sum(term[0].bracket(term[1]) for term in arg[0].terms())
        else:
            raise NotImplementedError

class SuperfunctionAlgebraUndirectedGraphOperation(SuperfunctionAlgebraOperation):
    """
    Operation on a SuperfunctionAlgebra defined by a UndirectedGraphVector.
    """
    def __init__(self, domain, codomain, graph_vector):
        """
        Initialize ``self``.
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
        terms = arg[0].terms()
        for e in graph.edges():
            new_terms = []
            for k in range(len(terms)):
                term0 = terms[k]
                left_sign = 1 if sum(term0[j].degree() for j in range(e[1])) % 2 == 0 else -1
                right_sign = 1 if sum(term0[j].degree() for j in range(e[0])) % 2 == 0 else -1
                for k in range(self._codomain.ngens()):
                    if 0 in term0:
                        continue
                    left_term = [f.copy() for f in term0]
                    left_term[e[0]] = left_term[e[0]].diff(evens[k])
                    left_term[e[1]] = left_sign * left_term[e[1]].diff(odds[k])
                    if not 0 in left_term:
                        new_terms.append(left_term)
                    right_term = [f.copy() for f in term0]
                    right_term[e[0]] = right_sign * right_term[e[0]].diff(odds[k])
                    right_term[e[1]] = right_term[e[1]].diff(evens[k])
                    if not 0 in right_term:
                        new_terms.append(right_term)
            terms = new_terms
        return sum(reduce(operator.mul, term) for term in terms)

    def __call__(self, *arg):
        """
        Return the evaluation of ``self`` at ``arg``.

        ASSUMPTION:

        Assumes that each factor in each term of ``arg`` is homogeneous.
        """
        result = 0
        for key in self._graph_vector._vector:
            g, sign = self._graph_vector.parent().basis().key_to_graph(key)
            c = self._graph_vector._vector[key]
            c *= sign
            result += c*self._act_with_graph(g, arg)
        return result

class SuperfunctionAlgebraSymmetricUndirectedGraphOperation(SuperfunctionAlgebraUndirectedGraphOperation):
    """
    Graded symmetric operation on a SuperfunctionAlgebra defined by a UndirectedGraphVector.
    """
    def __call__(self, *arg):
        """
        Return the evaluation of ``self`` at `arg``.

        ASSUMPTION:

        Assumes that each factor in each term of ``arg`` is homogeneous.
        """
        return super().__call__(arg[0].graded_symmetrization())
