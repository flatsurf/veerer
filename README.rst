Veerer
======

Veerer is a package for `SageMath <https://www.sagemath.org>`_ to deal with
veering triangulations of surfaces and their associated flat structures. It can
in particular be used to provide representatives of pseudo-Anosov mapping
classes of surfaces. The theoretical background is based on ideas of I. Agol
and F. Guéritaud, and is developed in:

M. Bell, V. Delecroix, V. Gadre, R. Gutiérrez-Romo, S. Schleimer,
"Coding Teichmüller flow using veering triangulations",
`arXiv:1909.00890 <https://arxiv.org/abs/1909.00890>`_.

Example
-------

As in `Flipper <https://github.com/MarkCBell/flipper>`_,
edges of a triangulation are labeled with the integers 0, 1, ..., n-1.
Each edge comes with an orientation and the edge opposite to ``i``
is labelled ``~i`` (so ``-1``, ``-2``, etc. are respectively opposite
edges to ``0``, ``1``, etc.). To input a triangulation, you must provide
a list of triangles, each triangle being a triple of oriented edges,
and a list of colours:: 

    >>> from veerer import *
    >>> T = VeeringTriangulation([(0, 1, 2), (-1, -2, -3)], [RED, RED, BLUE])
    >>> T.is_core()
    True

Since the above example is a core triangulation, it admits a flat realization.
One can be computed by taking the middle of the barycenter of the vertices of
the polytope parametrizing the flat structures::

    >>> F = T.flat_structure_middle()

If you use Veerer inside SageMath, the flat structure can be displayed with::

    >>> F.plot(vertical_train_track=True)
    >>> F.plot(horizontal_train_track=True)

Testing
-------

To run the SageMath doctests, install the package with pip, typically::

    $ sage -pip install -e .

and then run::

    $ sage -t --force-lib veerer/

Or::

    $ sage -t --force-lib veerer/my_file.py

Building the documentation
--------------------------

Go to the ``docs`` directory and then do::

    $ sage -sh
    $ make html

The documentation should be available under ``docs/build/`` as HTML pages.

Typically you might want to use ``veerer_demo.rst`` as a Jupyter notebook.
In order to convert ``veerer_demo.rst`` into ``veerer_demo.ipynb`` you need
to have available on your computer

- rst2latex python-docutils
- pdflatex 
- pandoc
- the Python package rst2ipynb
- the Python package nbconvert

Then do::

    $ export FILE_PREFIX="veerer_demo"
    $ rst2ipynb --kernel=sagemath veerer_demo.rst veerer_demo.ipynb

If you installed rst2ipynb using the ``--user`` option of pip, the executables
might be installed in ``$HOME/.local/bin``, in which case you should first make
the system aware of this via::

    $ PATH=$PATH:$HOME/.local/bin

Authors
-------

- Mark Bell
- Vincent Delecroix
- Saul Schleimer
