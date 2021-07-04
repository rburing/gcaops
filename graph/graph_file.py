from util.undirected_graph_sage import undirected_graph_from_encoding, undirected_graph_to_encoding
from util.directed_graph_sage import directed_graph_from_encoding, directed_graph_to_encoding
from abc import ABC, abstractmethod
import sqlite3
from math import ceil

class GraphFileView(ABC):
    def __init__(self, filename, num_vertices, num_edges):
        self._filename = filename
        self._num_vertices = num_vertices
        self._num_edges = num_edges
        self._con = sqlite3.connect(filename)
        cur = self._con.cursor()
        result = cur.execute('SELECT COUNT(name) FROM sqlite_master WHERE type = "table" AND name = "graphs"')
        if result.fetchone()[0] == 0:
            cur.execute('CREATE TABLE graphs (id INTEGER PRIMARY KEY, graph VARCHAR({0:d}))'.format(self._encoding_length()))
            self._con.commit()

    @abstractmethod
    def _encoding_to_graph(self, enc):
        pass

    @abstractmethod
    def _graph_to_encoding(self, g):
        pass

    @abstractmethod
    def _encoding_length(self):
        pass

    def __getitem__(self, index):
        cur = self._con.cursor()
        result = cur.execute('SELECT graph FROM graphs WHERE id = ?', (index + 1,)).fetchone()
        if result is not None:
            return self._encoding_to_graph(result[0])
        else:
            raise IndexError

    def index(self, g):
        encoding = self._graph_to_encoding(g)
        cur = self._con.cursor()
        result = cur.execute('SELECT id FROM graphs WHERE graph = ?', (encoding,)).fetchone()
        if result is not None:
            return result[0] - 1
        else:
            raise ValueError

    def __iter__(self):
        cur = self._con.cursor()
        for row in cur.execute('SELECT graph FROM graphs'):
            yield self._encoding_to_graph(row[0])

    def __len__(self):
        cur = self._con.cursor()
        result = cur.execute('SELECT COUNT(*) FROM graphs').fetchone()
        return result[0]

    def append(self, g):
        cur = self._con.cursor()
        cur.execute('INSERT INTO graphs (graph) VALUES (?)', (self._graph_to_encoding(g),))

    def commit(self):
        self._con.commit()

    # for pickling:

    def __getstate__(self):
        state_dict = self.__dict__.copy()
        del state_dict['_con']
        return state_dict

    def __setstate__(self, state_dict):
        self.__dict__.update(state_dict)
        self._con = sqlite3.connect(self._filename)

class UndirectedGraphFileView(GraphFileView):
    def _encoding_to_graph(self, enc):
        return undirected_graph_from_encoding(enc)

    def _graph_to_encoding(self, g):
        return undirected_graph_to_encoding(g)

    def _encoding_length(self):
        assert self._num_vertices <= 62
        num_bits = (self._num_vertices*(self._num_vertices-1)) // 2
        encoding_length = 1 + ceil(num_bits / 6.0) # graph6 length
        return encoding_length

class DirectedGraphFileView(GraphFileView):
    def _encoding_to_graph(self, enc):
        return directed_graph_from_encoding(enc)

    def _graph_to_encoding(self, g):
        return directed_graph_to_encoding(g)

    def _encoding_length(self):
        assert self._num_vertices <= 62
        num_bits = self._num_vertices * self._num_vertices
        encoding_length = len('&') + 1 + ceil(num_bits / 6.0) # digraph6 length
        return encoding_length

class UndirectedToDirectedGraphFileView:
    def __init__(self, filename):
        self._filename = filename
        self._con = sqlite3.connect(filename)
        cur = self._con.cursor()
        result = cur.execute('SELECT COUNT(name) FROM sqlite_master WHERE type = "table" AND name = "undirected_to_directed"')
        if result.fetchone()[0] == 0:
            cur.execute('CREATE TABLE undirected_to_directed (undirected_graph_id INTEGER, directed_graph_id INTEGER, coefficient INTEGER)')
            cur.execute('CREATE INDEX index_undirected ON undirected_to_directed(undirected_graph_id)')
            self._con.commit()

    def append(self, row):
        g_idx, h_idx, c = row
        cur = self._con.cursor()
        cur.execute('INSERT INTO undirected_to_directed (undirected_graph_id, directed_graph_id, coefficient) VALUES (?, ?, ?)', (1 + g_idx, 1 + h_idx, c))

    def commit(self):
        self._con.commit()
