# -*- coding: utf-8 -*-
import codecs
from setuptools import setup, Extension
from Cython.Build import cythonize

def readfile(filename):
    with codecs.open(filename,  encoding='utf-8') as f:
        return f.read()


extensions = [
    Extension('veerer.permutation', sources=['veerer/permutation.pyx'])
    ]

setup(
    name='veerer',
    version=readfile("VERSION"),
    description='A Python module to manipulate Veering triangulations and their associated flat structures',
    long_description=readfile("README.rst"),
    author='Mark Bell, Vincent Delecroix and Saul Schleimer',
    author_email='vincent.delecroix@u-bordeaux.fr',
    url='https://github.com/flatsurf/veerer',
    packages=['veerer/', 'veerer/polyhedron'],
    ext_modules=cythonize(extensions),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Mathematics',
        ]
)
