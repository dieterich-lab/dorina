# vim: set fileencoding=utf-8 :

import unittest
from os import path
from argparse import Namespace
from dorina import config
from dorina import utils

datadir = path.join(path.dirname(path.abspath(__file__)), 'data')

class TestListDataWithoutOptions(unittest.TestCase):
    def test_get_genomes(self):
        """Test utils.get_genomes()"""
        expected = {
            'mammals': {
                'h_sapiens': {
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

        got = utils.get_genomes(datadir=datadir)
        self.assertEqual(expected, got)


    def test_get_regulators(self):
        """Test utils.get_regulators()"""
        basedir = path.join(datadir, 'regulators', 'mammals', 'h_sapiens', 'hg19')
        scifi = utils.parse_experiment(path.join(datadir, basedir, 'RBP', 'scifi.json'))
        fake01 = utils.parse_experiment(path.join(datadir, basedir, 'miRNA', 'fake01.json'))
        fake02 = utils.parse_experiment(path.join(datadir, basedir, 'miRNA', 'fake02.json'))
        expected = {
            'mammals': {
                'h_sapiens': {
                    'hg19': {
                        'RBP': {
                            'scifi': scifi
                        },
                        'miRNA': {
                            'fake01': fake01,
                            'fake02': fake02
                        },
                    }
                }
            }
        }

        got = utils.get_regulators(datadir=datadir)
        self.assertEqual(expected, got)


    def test_parse_experiment(self):
        """Test utils.parse_experiment()"""
        expected = {
            'summary': 'Experimental summary',
            'description': 'Long description here',
            'methods': 'Experimental methods section',
            'credits': 'Credits',
            'references': [
                { 'title': 'A very important publication',
                  'authors': ['Jules Verne', 'Orson Wells'],
                  'pages': '23-42', 'journal': 'Annals of Science Fiction',
                  'year': '1870',
                  'pubmed': 'http://www.ncbi.nlm.nih.gov/pubmed/12345678'
                }
            ]
        }

        basedir = path.join(datadir, 'regulators', 'mammals', 'h_sapiens', 'hg19')
        got = utils.parse_experiment(path.join(datadir, basedir, 'RBP', 'scifi.json'))
        self.assertEqual(expected, got)


    def test_get_genome_by_name(self):
        """Test utils.get_genome_by_name()"""
        got = utils.get_genome_by_name("invalid", datadir=datadir)
        self.assertIsNone(got)

        expected = path.join(datadir, 'genomes', 'mammals', 'h_sapiens', 'hg19')
        got = utils.get_genome_by_name("hg19", datadir=datadir)
        self.assertEqual(expected, got)


    def test_get_regulator_by_name(self):
        """Test utils.get_regulator_by_name()"""
        got = utils.get_regulator_by_name("invalid", datadir=datadir)
        self.assertIsNone(got)

        expected = path.join(datadir, 'regulators', 'mammals', 'h_sapiens', 'hg19', 'RBP', 'scifi')
        got = utils.get_regulator_by_name("scifi", datadir=datadir)
        self.assertEqual(expected, got)
