from abc import ABC, abstractmethod

class GraphVector(ABC):
    """
    Vector representing a linear combination of graphs.
    """
    @abstractmethod
    def __repr__(self):
        """
        Return a string representation of this graph vector.
        """
        pass

    @abstractmethod
    def parent(self):
        """
        Return the parent GraphModule that this graph vector belongs to.
        """
        pass

    @abstractmethod
    def copy(self):
        """
        Return a copy of this graph vector.
        """
        pass

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

    @abstractmethod
    def __pos__(self):
        """
        Return a copy of this graph vector.
        """
        pass

    @abstractmethod
    def __neg__(self):
        """
        Return the negative of this graph vector.
        """
        pass

    @abstractmethod
    def __add__(self, other):
        """
        Return this graph vector added to ``other``.
        """
        pass

    @abstractmethod
    def __radd__(self, other):
        """
        Return ``other`` added to this graph vector.
        """
        pass

    @abstractmethod
    def __mul__(self, other):
        """
        Return this graph vector multiplied by ``other``.
        """
        pass

    @abstractmethod
    def __rmul__(self, other):
        """
        Return ``other`` multiplied by this graph vector.
        """
        pass

    @abstractmethod
    def __eq__(self, other):
        """
        Return ``True`` if this graph vector is equal to ``other`` and ``False`` otherwise.
        """
        pass

    @abstractmethod
    def bi_gradings(self):
        """
        Return the set of tuples ``(v,e)`` such that this graph vector contains terms with ``v`` vertices and ``e`` edges.
        """
        pass

    @abstractmethod
    def insertion(self, position, other):
        """
        Return the insertion of ``other`` into this graph vector at the vertex ``position``.
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

class GraphModule(ABC):
    """
    Module spanned by graphs.
    """
    @abstractmethod
    def base_ring(self):
        """
        Return the base ring of this module.
        """
        pass

    @abstractmethod
    def basis(self):
        """
        Return the basis of this module.
        """
        pass

    @abstractmethod
    def __repr__(self):
        """
        Return a string representation of this module.
        """
        pass

    @abstractmethod
    def zero(self):
        """
        Return the zero vector in this module.
        """
        pass

    @abstractmethod
    def __call__(self, arg):
        """
        Convert ``arg`` into an element of this module.
        """
        pass
