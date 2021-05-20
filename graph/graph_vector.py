from abc import ABC, abstractmethod

class GraphVector(ABC):
    """
    Vector representing a linear combination of graphs.
    """
    @abstractmethod
    def __iter__(self):
        """
        Returns an iterator over this graph vector, yielding tuples of the form ``(coeff, graph)``.
        """
        pass
