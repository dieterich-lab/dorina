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
utils = dorina.utils.DorinaUtils(datadir)

class TestListDataWithoutOptions(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None

    def test_genome_path_by_name(self):
        """Test Genome.path_by_name()"""
        got = Genome.path_by_name("invalid")
        self.assertIsNone(got)

        expected = path.join(datadir, 'genomes', 'h_sapiens', 'hg19')
        got = Genome.path_by_name(utils.genomes, "hg19")
        self.assertEqual(expected, got)

    def test_get_genes(self):
        """Test utils.get_genes()"""
        expected = ['gene01.01', 'gene01.02']
        got = utils.get_genes('hg19')
        self.assertEqual(expected, got)

        expected = []
        got = utils.get_genes('invalid')
        self.assertEqual(expected, got)
