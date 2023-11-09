# ********************************************************************
#  This file is part of veerer
#
#        Copyright (C) 2023 Julian Rüth
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ********************************************************************

try:
  input("Are you sure you are on the master branch which is identical to origin/master? [ENTER]")
except KeyboardInterrupt:
  sys.exit(1)

$PROJECT = 'veerer'

$ACTIVITIES = [
    'version_bump',
    'tag',
    'push_tag',
    'pypi',
    'ghrelease',
]

$RELEASE_YEAR = $RELEASE_DATE.year

$VERSION_BUMP_PATTERNS = [
    ('veerer/version.py', r"version =", "version = \"$VERSION\""),
    ('setup.py', r"    version=", "    version=\"$VERSION\","),
    ('doc/source/conf.py', r'copyright = ', "copyright = \"2016-$RELEASE_YEAR, Mark Bell, Vincent Delecroix, Saul Schleimer\""),
]

$PUSH_TAG_REMOTE = 'git@github.com:flatsurf/veerer.git'

$PYPI_BUILD_COMMANDS = ['sdist']
$PYPI_NAME = "veerer"

$GITHUB_ORG = 'flatsurf'
$GITHUB_REPO = 'veerer'
