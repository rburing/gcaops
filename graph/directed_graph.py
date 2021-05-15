from collections.abc import MutableSequence
from util.permutation import selection_sort

class DirectedGraph:
    """
    Directed graph with vertices labeled by natural numbers and an ordered set of edges.
    """
    def __init__(self, num_vertices, edges):
        """
        Initialize this directed graph.

        INPUT:

        - ``num_vertices`` -- a natural number, the number of vertices

        - ``edges`` -- a list of tuples of natural numbers
        """
        if not num_vertices >= 0:
            raise ValueError('num_vertices must be a natural number')
        self._num_vertices = num_vertices
        if not isinstance(edges, MutableSequence) or not all(isinstance(edge, tuple) for edge in edges):
            raise ValueError('Format of edges {} not recognized'.format(edges))
        for (source,target) in edges:
            if source >= num_vertices or target >= num_vertices:
                raise ValueError('Vertex labels must be between 0 and the number of vertices')
        self._edges = edges

    def __repr__(self):
        """
        Return a string representation of this graph.
        """
        return 'DirectedGraph({}, {})'.format(self._num_vertices, self._edges)

    def __len__(self):
        """
        Return the number of vertices of this graph.
        """
        return self._num_vertices

    def __eq__(self, other):
        """
        Return ``True`` if this graph equals ``other``.
        """
        return isinstance(other, self.__class__) and self._num_vertices == other._num_vertices and self._edges == other._edges

    def edges(self):
        """
        Return the list of edges of this graph.
        """
        return self._edges

    def canonicalize_edges(self):
        """
        Lexicographically order the edges of this graph and return the sign of that edge permutation.
        """
        return selection_sort(self._edges)

    def relabeled(self, relabeling):
        """
        Return a vertex relabeling of this graph.
        """
        new_edges = [(relabeling[a], relabeling[b]) for (a,b) in self._edges]
        return __class__(self._num_vertices, new_edges)
