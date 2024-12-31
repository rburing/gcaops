# `gcaops`: Graph Complex Action on Poisson Structures

Python package implementing the action of Kontsevich's graph complex(es) on Poisson structures.
This package is designed to be used in [SageMath](https://www.sagemath.org/) version 9.5 or later.

## Installation in SageMath

1. Download the `gcaops` source code as a ZIP file and extract it to a directory such as `/path/to/gcaops-master`.
2. In a terminal (e.g. the SageMath Shell on Windows), run the following:
    ```
    $ sage -pip install --upgrade /path/to/gcaops-master/
    ```
    This completes the installation.
3. It is optional but highly recommended to configure a default directory where data (such as lists of graphs) can be stored, so it doesn't have to be re-computed each time.

    This can be done by setting the environment variable `GCAOPS_DATA_DIR` to the path you desire, before starting SageMath.
    A convenient way to achieve this is by adding a line such as the following to SageMath's [shell script](https://doc.sagemath.org/html/en/reference/repl/startup.html#the-sagerc-shell-script) `sagerc`:
    ```
    export GCAOPS_DATA_DIR='/home/sage/Documents/gcaops_data/'
    ```
    Be warned that this directory can grow large. If no directory is configured, then graphs are only stored in memory (which may be limiting).
4. It is optional but convenient to enable the importing of all names from the `gcaops` package (e.g. `UndirectedGraphComplex`) into the global namespace of every SageMath session, so that the functionality can be used immediately.

    This can be done by adding the following line to SageMath's [startup script](https://doc.sagemath.org/html/en/reference/repl/startup.html#the-init-sage-script) `init.sage`:
    ```
    from gcaops.all import *
    ```

## Usage in SageMath

In a SageMath session, import the package and use it. For instance:

```
sage: from gcaops.graph.undirected_graph_complex import UndirectedGraphComplex
sage: GC = UndirectedGraphComplex(QQ, implementation='vector', sparse=True)
sage: GC.cohomology_basis(4,6)
[1*UndirectedGraph(4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])]
```

Extensive examples of the use of this software are contained in the author's PhD dissertation: [*The action of Kontsevich's graph complex on Poisson structures and star products: an implementation*](https://doi.org/10.25358/openscience-9274).

## Installation in a virtual Python environment (no SageMath installation required)

Create and activate a virtual environment:

    cd /path/to/gcaops-master/
    python3 -m venv venv_gcaops
    . venv_gcaops/bin/activate

Install the package in the virtual environment:

    pip install ".[passagemath]"

This automatically installs the modularized parts of the Sage library that are
needed by the package. (These modularized distributions are provided by
https://github.com/passagemath.)
