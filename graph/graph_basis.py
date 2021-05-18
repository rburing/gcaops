from abc import ABC, abstractmethod

class GraphBasis(ABC):
    """
    Basis of a module spanned by graphs.

    A basis consists of keys ``(v,e,index,...)`` where ``(v,e,index)`` identifies the isomorphism class of the graph.
    """
    graph_class = None.__class__

    @abstractmethod
    def _generate_graphs(self, bi_grading):
        """
        Return a list of all the graphs in this basis in the given ``bi_grading``.
        """
        pass

    @abstractmethod
    def graphs(self, vertices, edges):
        """
        Return the list of graphs in this basis with the given amount of ``vertices`` and ``edges``.
        """
        pass

    @abstractmethod
    def cardinality(self, vertices, edges):
        """
        Return the number of graphs in this basis with the given amount of ``vertices`` and ``edges``.
        """
        pass

    @abstractmethod
    def graph_to_key(self, graph):
        """
        Return a tuple consisting of the key in this basis and the sign factor such that ``graph`` equals the sign times the graph identified by the key.

        INPUT:

        - ``graph`` -- a graph
        """
        pass

    @abstractmethod
    def key_to_graph(self, key):
        """
        Return a tuple consisting of a graph and the sign factor such that the sign times the graph equals the graph identified by the key.

        INPUT:

        - ``key`` -- a key in this basis
        """
        pass

    def __repr__(self):
        """
        Return a string representation of this basis.
        """
        return 'Basis consisting of graphs'

