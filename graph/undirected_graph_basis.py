class UndirectedGraphBasis:
    """
    Basis of a module of undirected graphs.
    """
    def __init__(self, graph_to_index, index_to_graph):
        """
        Initialize ``self``.

        INPUT:

        - ``graph_to_index`` -- a function that takes an ``UndirectedGraph`` and returns the index in the basis and a sign

        - ``index_to_graph`` -- a function that takes an index in the basis and returns an ``UndirectedGraph`` and a sign
        """
        self.graph_to_index = graph_to_index
        self.index_to_graph = index_to_graph

    def __repr__(self):
        return 'Basis of undirected graphs'
