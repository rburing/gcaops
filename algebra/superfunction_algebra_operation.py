from abc import ABC, abstractmethod
from .tensor_product import TensorProduct

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
        return 'Operation on {} of arity {}'.format(self._codomain, self._domain.nfactors())

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
