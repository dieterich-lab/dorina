#!/usr/bin/env python
# -*- coding: utf-8

from __future__ import unicode_literals
import unittest
from os import path

from dorina.genome import Genome


class TestListDataWithoutOptions(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.datadir = path.join(path.dirname(path.abspath(__file__)), 'data')
        self.Genome = Genome.init(self.datadir)

    def TearDown(self):
        self.maxDiff = None
        self.datadir = None
        self.Genome = None

    def test_genomes(self):
        """Test Genome.all"""
        expected = {
            'h_sapiens': {
                'id': 'h_sapiens',
                'label': 'Human',
                'scientific': 'Homo sapiens',
                'weight': 10,
                'assemblies': {
                    'hg19': {
                        'all': True,
                        'cds': True,
                        '3_utr': True,
                        '5_utr': True,
                        'intron': True,
                        'intergenic': True
                    }
                }
            }
        }

        got = Genome.all()
        self.assertEqual(expected, got)

    def test_genome_path_by_name(self):
        """Test Genome.path_by_name()"""
        with self.assertRaises(ValueError):
            got = Genome.path_by_name("invalid")

        expected = path.join(self.datadir, 'genomes', 'h_sapiens', 'hg19')
        got = Genome.path_by_name("hg19")
        self.assertEqual(expected, got)

    def test_get_genes(self):
        """Test Genome.get_genes()"""
        expected = ['gene01.01', 'gene01.02']
        got = Genome.get_genes('hg19')
        self.assertEqual(expected, got)
