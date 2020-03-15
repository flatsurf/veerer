#!/usr/bin/env python

import subprocess
import sys
import os.path

directory = os.path.dirname(__file__)

import veerer.env

tests = ["permutation.py", "triangulation_comparison.py", "triangulation_isomorphism.py", "flip.py"]

if veerer.env.ppl is not None:
    tests.extend(["geometric_polytope.py"])

for test in tests:
    test = os.path.join(directory, test)
    subprocess.run([sys.executable, test] + sys.argv[1:])
