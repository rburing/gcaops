from collections.abc import MutableSequence
from util.permutation import selection_sort

class UndirectedGraph:
    """
    Undirected graph with vertices labeled by natural numbers and an ordered set of edges.
    """
    def __init__(self, num_vertices, edges):
        """
        Initialize ``self``.

        INPUT:

        - ``num_vertices`` -- a natural number, the number of vertices

        - ``edges`` -- a list of tuples of natural numbers
        """
        if not num_vertices >= 0:
            raise ValueError('num_vertices must be a natural number')
        self._num_vertices = num_vertices
        if not isinstance(edges, MutableSequence) or not all(isinstance(edge, tuple) for edge in edges):
            raise ValueError('Format of edges {} not recognized'.format(edges))
        self._edges = edges
        # canonicalize representation of individual edges
        for k in range(len(self._edges)):
            first, last = self._edges[k]
            if first >= num_vertices or last >= num_vertices:
                raise ValueError('Vertex labels must be between 0 and the number of vertices')
            if first > last:
                self._edges[k] = (last, first)

    def __repr__(self):
        """
        Return a string representation of ``self``.
        """
        return 'UndirectedGraph({}, {})'.format(self._num_vertices, self._edges)

    def edges(self):
        """
        Return the list of edges of ``self``.
        """
        return self._edges

    def canonicalize_edges(self):
        """
        Lexicographically order the edges of ``self`` and return the sign of that edge permutation.
        """
        return selection_sort(self._edges)
