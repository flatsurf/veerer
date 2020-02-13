Veerer
======

veerer is a Python module to deal with veering triangulations of surfaces and
their associated flat structures. It can in particular be used to provides
representatives of pseudo-Anosov mapping classes of surfaces. The theoretical
background is developed in M. Bell, V. Delecroix, V. Gadre, R. Gutiérrez-Romo,
S. Schleimer "Coding Teichmüller flow using veering triangulations"
(`<https://arxiv.org/abs/1909.00890>`_) and is based on ideas of
I. Agol and F. Guéritaud.

To install the module you need Python (preferably version 3 but Python 2 is
supported). Computations involving polytopes are only available if the Python
module pplpy is available (see https://gitlab.com/videlec/pplpy). It is
installed by default with SageMath from version 8.7. Additional features are
available if this module is used inside SageMath (http://www.sagemath.org/).

Example
-------

As with `<https://github.com/MarkCBell/flipper>`_, edges of a triangulation are
labeled with the integers 0, 1, ..., n-1. Each edge come with an orientation and
the edge opposite to i is labelled ~i (that is -1 is the opposite of 0, -2 of 1,
etc). To input a triangulation, you must provide a list of triangles, each
triangle beeing a triple of oriented edges, and a list of colours::

    >>> from veerer import *
    >>> T = VeeringTriangulation([(0,1,2), (-1,-2,-3)], [RED, RED, BLUE])
    >>> T.is_core()
    True

Since the above example is a core triangulation it admits flat realization. One
can be computed by taking the middle of the barycenter of the vertices of
the polytope parametrizing the flat structures::

    >>> F = T.flat_structure_middle()

If you veerer inside SageMath, the flat structure could be displayed with::

    >>> F.plot(vertical_train_track=True)
    >>> F.plot(horizontal_train_track=True)

Testing
-------

To run the SageMath doctests, install the module with pip, typically::

    $ sage -pip install . --user --force-reinstall

and then run::

    $ sage -t --force-lib veerer/

Or::

    $ sage -t --force-lib veerer/my_file.py

Building documentation
----------------------

Go to the docs repository and then do::

    $ sage -sh
    $ make html

The documentation should be available under docs/build/ as HTML pages.

Typically you might want to use veerer_demo.rst as a Jupyter notebook. In
order to do the conversion you need to have available on your computer

- rst2latex python-docutils
- pdflatex 
- pandoc
- the Python module rst2ipynb
- the Python module nbconvert

Then do::

    $ export FILE_PREFIX="veerer_demo"
    $ rst2ipynb --kernel=sagemath veerer_demo.rst veerer_demo.ipynb

If you did install rst2ipynb using the `--user` option of pip the executables
are possibly installed in $HOME/.local/bin. In which case you should first make
the system aware of this via::

    $ PATH=$PATH:$HOME/.local/bin

Authors
-------

- Mark Bell
- Vincent Delecroix
- Saul Schleimer
