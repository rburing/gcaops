from .graph_vector import GraphVector, GraphModule
from abc import ABC, abstractmethod

class GraphCochain(GraphVector):
    """
    Cochain of a GraphComplex.
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
