from util.misc import keydefaultdict
from util.permutation import selection_sort
from util.nauty import nauty_canonicalize, nauty_generate, nauty_has_odd_automorphism

class UndirectedGraphBasis:
    """
    Basis of a module of undirected graphs.
    """
    def __init__(self):
        """
        Initialize ``self``.
        """
        self._graphs = keydefaultdict(lambda key: list(filter(lambda g: not nauty_has_odd_automorphism(g), nauty_generate(*key))))

    def graph_to_index(self, graph):
        """
        Return the index in ``self`` and the sign corresponding to the input ``graph``.

        INPUT:

        - ``graph`` - an UndirectedGraph
        """
        g, sign = nauty_canonicalize(graph)
        v, e = len(g), len(g.edges())
        try:
            index = self._graphs[v,e].index(g)
            return (v,e,index), sign
        except ValueError:
            return None, 1

    def index_to_graph(self, index):
        """
        Return the UndirectedGraph and the sign corresponding to the index ``index`` in ``self``.

        INPUT:

        - ``index`` - an index in ``self``
        """
        v, e, idx = index
        return self._graphs[v,e][idx], 1

    def __repr__(self):
        """
        Return a string representation of ``self``.
        """
        return 'Basis of undirected graphs'
