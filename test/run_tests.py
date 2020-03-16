#!/usr/bin/env python

import subprocess
import sys
import os.path

directory = os.path.dirname(__file__)

import veerer.env

tests = ["permutation.py", "triangulation_comparison.py", "triangulation_reconstruction.py", "triangulation_relabel.py", "triangulation_isomorphism.py", "flip_sequence_arithmetic.py"]

if veerer.env.ppl is not None:
    tests.extend(["flip.py", "geometric_polytope.py"])

output = 0
for test in tests:
    test = os.path.join(directory, test)
    output += subprocess.call([sys.executable, test] + sys.argv[1:])

sys.exit(output)
