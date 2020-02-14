# -*- coding: utf-8 -*-
import os
import codecs
from distutils.core import setup
from distutils.cmd import Command

def readfile(filename):
    with codecs.open(filename,  encoding='utf-8') as f:
        return f.read()

class TestCommand(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import subprocess

        if subprocess.call(['sage', '-t', '--force-lib', 'veerer/']):
            raise SystemExit("Doctest failures")


setup(
    name='veerer',
    version=readfile("VERSION"),
    description='A Python module to manipulate Veering triangulations and their associated flat structures',
    long_description=readfile("README.rst"),
    author='Mark Bell, Vincent Delecroix and Saul Schleimer',
    author_email='vincent.delecroix@u-bordeaux.fr',
    url='https://framagit.org/saulsch/Veerer',
    # Remember to update these if the directory structure changes.
    packages=['veerer/'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Mathematics',
        ],
    cmdclass = {'test': TestCommand}
)
