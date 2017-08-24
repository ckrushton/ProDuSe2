#! /usr/bin/env python

from setuptools import setup
import re

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
	author="Nolan Martin",
	author_email="ncm3@sfu.ca",
	maintainer="Christopher Rushton",
	maintainer_email="ckrushto@sfu.ca",
	packages=["ProDuSe"],
	install_requires=[
		"pysam",
		"sortedcontainers",
		"networkx",
		"ruamel.yaml",
		"jmespath",
		"asciimatics"],
	license="GNU GPLv3",
	scripts=["bin/produse"]
	)
