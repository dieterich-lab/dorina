# vim: set fileencoding=utf-8 :

import unittest
from os import path
from argparse import Namespace
import json

import dorina
from dorina import config
from dorina.genome import Genome
from pybedtools import BedTool

datadir = path.join(path.dirname(path.abspath(__file__)), 'data')

class TestListDataWithoutOptions(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None

    def test_get_genes(self):
        """Test Genome.get_genes()"""
        expected = ['gene01.01', 'gene01.02']
        got = Genome.get_genes('hg19')
        self.assertEqual(expected, got)
