# ****************************************************************************
#  This file is part of veerer
#
#       Copyright (C) 2019-2023 Vincent Delecroix
#                          2023 Julian Rüth
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# ****************************************************************************
from setuptools import setup, Extension
from Cython.Build import cythonize


with open("README.rst") as f:
    long_description = f.read()


extensions = [Extension("veerer.permutation", sources=["veerer/permutation.pyx"])]

setup(
    name="veerer",
    version="0.2.0",
    description="A Python module to manipulate Veering triangulations and their associated flat structures",
    long_description=long_description,
    author="Mark Bell, Vincent Delecroix, and Saul Schleimer",
    author_email="contact@flatsurf.org",
    url="https://github.com/flatsurf/veerer",
    packages=["veerer", "veerer.polyhedron"],
    package_data={"veerer": ["permutation.pyx"]},
    license = 'GPL 2.0+',
    ext_modules=cythonize(extensions),
    install_requires=[],
    setup_requires=["Cython"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
)
