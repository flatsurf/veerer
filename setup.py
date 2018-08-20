import os
from distutils.core import setup

with open(os.path.join(os.path.dirname(__file__), 'README')) as f:
    README = f.read()

setup(
    name='veerer',
    version='0.0.1.beta',
    description='A Python module to manipulate Veering triangulations and their associated flat structures',
    long_description=README,
    author='Mark Bell, Vincent Delecroix and Saul Schleimer',
    author_email='vincent.delecroix@u-bordeaux.fr',
    url='https://framagit.org/saulsch/Veerer',
    # Remember to update these if the directory structure changes.
    packages=['veerer/'],
#    install_requires=[
#        'pplpy',
#        'surface_dynamics'
#        ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Mathematics',
        ],
)

