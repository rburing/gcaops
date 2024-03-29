r"""
Graph complex
"""
from abc import abstractmethod
from .graph_vector import GraphVector, GraphModule

class GraphCochain(GraphVector):
    """
    Cochain of a :class:`GraphComplex`.
    """
    @abstractmethod
    def bracket(self, other):
        """
        Return the graph Lie bracket of this graph cochain with ``other``.
        """
        pass

    @abstractmethod
    def differential(self):
        """
        Return the graph differential of this graph cochain.
        """
        pass

class GraphComplex(GraphModule):
    """
    Graph complex.
    """
    pass
