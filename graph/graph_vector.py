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

    @abstractmethod
    def __len__(self):
        """
        Return the number of graphs with nonzero coefficients in this graph vector.
        """
        pass

    def plot(self, **options):
        """
        Return a plot of this graph vector.
        """
        from sage.misc.latex import latex
        from sage.plot.plot import graphics_array
        from sage.plot.text import text
        my_options = {'xmin': 0.0, 'xmax': 3.0, 'ymin': -1.0, 'ymax': 1.0, 'figsize': [1.0, 1.0], 'aspect_ratio': 1.0, 'axes': False}
        my_options.update(options)
        fontsize = my_options.pop('fontsize') if 'fontsize' in my_options else 16
        ncols = my_options.pop('ncols') if 'ncols' in my_options else 3
        label = lambda c: text('${}' + (latex(c) if str(c)[0] == '-' else '+' + latex(c)) + '$', (0.4,0.0), fontsize=fontsize, axes=False)
        return graphics_array([g.plot(**my_options) + label(c) for (c,g) in self], ncols=ncols)

    def show(self, **options):
        """
        Show this graph.
        """
        my_options = {'figsize': [3 * 3.5, 2.0 + len(self)//3 * 2.0], 'aspect_ratio': 1.0}
        my_options.update(options)
        from sage.graphs.graph_plot import graphplot_options
        plot_options = {k: my_options.pop(k) for k in graphplot_options if k in my_options}
        return self.plot(**plot_options).show(**my_options)
