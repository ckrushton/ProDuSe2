#! /usr/bin/env python
from setuptools import setup

# Imports version number
try:
    from ProDuSe.__version import __version__ as version
except ImportError:
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
