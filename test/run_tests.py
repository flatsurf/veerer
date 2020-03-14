#!/usr/bin/env python

import subprocess
import sys

import veerer.env

tests = ["permutation.py", "triangulation_isomorphism.py", "flip.py"]

if veerer.env.ppl is not None:
    tests.extend(["geometric_polytope.py"])

for test in tests:
    subprocess.run([sys.executable, test] + sys.argv[1:])
