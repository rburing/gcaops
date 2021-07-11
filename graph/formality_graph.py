from collections.abc import MutableSequence
from util.permutation import selection_sort

class FormalityGraph:
    """
    Directed graph with an ordered set of edges, and vertices labeled by natural numbers, the first of which are ordered ground vertices without outgoing edges.
    """
    def __init__(self, num_ground_vertices, num_aerial_vertices, edges):
        """
        Initialize this formality graph.

        INPUT:

        - ``num_ground_vertices`` -- a natural number, the number of ground vertices

        - ``num_aerial_vertices`` -- a natural number, the number of aerial vertices

        - ``edges`` -- a list of tuples of natural numbers
        """
        if not num_ground_vertices >= 0:
            raise ValueError('num_ground_vertices must be a natural number')
        self._num_ground_vertices = num_ground_vertices
        if not num_aerial_vertices >= 0:
            raise ValueError('num_aerial_vertices must be a natural number')
        self._num_aerial_vertices = num_aerial_vertices
        if not isinstance(edges, MutableSequence) or not all(isinstance(edge, tuple) for edge in edges):
            raise ValueError('Format of edges {} not recognized'.format(edges))
        num_vertices = num_ground_vertices + num_aerial_vertices
        for (source,target) in edges:
            if source >= num_vertices or target >= num_vertices:
                raise ValueError('Vertex labels must be between 0 and the total number of vertices')
            if source < num_ground_vertices:
                raise ValueError('Ground vertices must not have any outgoing edges')
        self._edges = edges
        self._vertex_positions = None

    def __repr__(self):
        """
        Return a string representation of this graph.
        """
        return 'FormalityGraph({}, {}, {})'.format(self._num_ground_vertices, self._num_aerial_vertices, self._edges)

    def __len__(self):
        """
        Return the number of vertices of this graph.
        """
        return self._num_aerial_vertices + self._num_ground_vertices

    def __eq__(self, other):
        """
        Return ``True`` if this graph equals ``other``.
        """
        return isinstance(other, self.__class__) and self._num_ground_vertices == other._num_ground_vertices and self._num_aerial_vertices == other._num_aerial_vertices and self._edges == other._edges

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
        if any((v < self._num_ground_vertices) != (w < self._num_ground_vertices) for (v,w) in relabeling.items()):
            raise ValueError('Relabeling must map aerial vertices to aerial vertices')
        new_edges = [(relabeling[a], relabeling[b]) for (a,b) in self._edges]
        return __class__(self._num_ground_vertices, self._num_aerial_vertices, new_edges)

    def out_degrees(self):
        """
        Return the tuple of out-degrees of vertices of this graph.
        """
        degrees = [0 for i in range(self._num_ground_vertices + self._num_aerial_vertices)]
        for (a,b) in self._edges:
            degrees[a] += 1
        return tuple(degrees)

    def get_pos(self):
        """
        Return the dictionary of positions of vertices in this graph (used for plotting).
        """
        return self._vertex_positions

    def set_pos(self, new_pos):
        """
        Set the positions of vertices in this graph (used for plotting).
        """
        self._vertex_positions = new_pos

    def plot(self, **options):
        """
        Return a plot of this graph.
        """
        from sage.graphs.digraph import DiGraph
        from sage.graphs.graph_plot import GraphPlot
        g = DiGraph([(a,b,i) for (i,(a,b)) in enumerate(self.edges())])
        ground_pos = { v : (v,0) for v in range(self._num_ground_vertices)}
        aerial_vertices = range(self._num_ground_vertices, self._num_ground_vertices + self._num_aerial_vertices)
        vertex_positions = self.get_pos()
        if vertex_positions:
            g.set_pos(vertex_positions)
            plot = GraphPlot(graph=g, options=options).plot()
            pos = g.get_pos()
        else:
            # NOTE: naively retries plotting until all aerial vertices are above the real axis
            while True:
                new_pos = {}
                new_pos.update(ground_pos)
                g.set_pos(new_pos)
                plot = GraphPlot(graph=g, options=options).plot()
                pos = g.get_pos()
                if all(pos[v][1] > 0 for v in aerial_vertices):
                    break
        if options.get('save_pos', False):
            self.set_pos(pos)
        return plot

    def show(self, **options):
        """
        Show this graph.
        """
        from sage.graphs.graph_plot import graphplot_options
        plot_options = {k: options.pop(k) for k in graphplot_options if k in options}
        return self.plot(**plot_options).show(**options)
