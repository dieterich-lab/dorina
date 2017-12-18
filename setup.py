#!/usr/bin/env python
# -*- coding: utf-8
import os
from setuptools import setup, find_packages
from dorina import __version__


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="dorina",
    version=__version__,
    author="Kai Blin",
    author_email="kai.blin@age.mpg.de",
    maintainer="Thiago Britto Borges",
    maintainer_email="thiago.brittoborges@uni-heidelberg.de",
    description="database of posttranscriptional regulatory elements",
    keywords="bioinformatics",
    url="https://github.com/dieterich-lab/dorina",
    packages=find_packages(),
    package_data={'dorina.config': ['*.cfg']},
    include_package_data=True,
    install_requires=['Cython', 'pybedtools'],
    tests_require=['nose'],
    license=read('LICENSE'),
    long_description=read('README.md'),

    entry_points='''
    [console_scripts]
    dorina=dorina.__main__.cli''',
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3.6",

    ],
)
