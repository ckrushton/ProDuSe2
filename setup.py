#! /usr/bin/env python
from setuptools import setup
import re
import sys

# Check the python version, and ensure it is supported
pythonVer = ".".join(list(str(x) for x in sys.version_info[0:2]))
if float(pythonVer) < 3.4:
    sys.stderr.write("ERROR: Python 3.4 or newer is required to install and run ProDuSe\n")
    exit(1)

# Imports version number
VERSIONFILE = "ProDuSe/__version.py"
verstrline = open(VERSIONFILE, "rt").read()
verRegex = r"^__version__ = ['\"]([^'\"]*)['\"]"
currentVer = re.search(verRegex, verstrline, re.M)
if currentVer:
    version = currentVer.group(1)
else:
    version = "Unknown"

setup(
    name="ProDuSe",
    version=version,
    description="Variant caller for semi-degenerate barcoded adapter libraries",
    author="Nolan",
    author_email="innovate.invent@gmail.com",
    maintainer="Christopher Rushton",
    maintainer_email="ckrushto@sfu.ca",
    packages=["ProDuSe"],
    python_requires=">3.4",
    install_requires=[
        "pysam",
        "sortedcontainers",
        "networkx",
        "configutator",
        "clipoverlap",
        "CigarIterator"],
    license="GNU GPLv3",
    scripts=["bin/produse"]
    )
