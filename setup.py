#!/usr/bin/env python
# -*- coding: utf-8
import os
from setuptools import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


version = 'undefined'
for line in open(os.path.join('dorina', '__init__.py'), 'r'):
    if line.startswith('version'):
        exec (line.strip())

setup(
    name="dorina",
    version=version,
    author="Kai Blin",
    author_email="kai.blin@age.mpg.de",
    description=("database of posttranscriptional regulatory elements"),
    keywords="bioinformatics",
    url="https://bioinf-redmine.age.mpg.de/projects/dorina-2",
    packages=['dorina', 'dorina.config'],
    install_requires=['Cython>=0.20.1', 'pybedtools>=0.6.4'],
    tests_require=['nose'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
