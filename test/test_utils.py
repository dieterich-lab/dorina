# vim: set fileencoding=utf-8 :

import unittest
from os import path
from argparse import Namespace
from dorina import config
from dorina import utils

datadir = path.join(path.dirname(path.abspath(__file__)), 'data')

class TestListData(unittest.TestCase):
    def setUp(self):
        self.old_config = config.get_config()
        options = Namespace()
        options.data = Namespace()
        options.data.path = datadir
        config.set_config(options)

    def tearDown(self):
        config.set_config(self.old_config)

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

        got = utils.get_genomes()
        self.assertEqual(expected, got)


    def test_get_regulators(self):
        """Test utils.get_regulators()"""
        scifi = utils.parse_experiment(path.join(datadir, 'regulators', 'RBP', 'scifi.json'))
        fake01 = utils.parse_experiment(path.join(datadir, 'regulators', 'miRNA', 'fake01.json'))
        fake02 = utils.parse_experiment(path.join(datadir, 'regulators', 'miRNA', 'fake02.json'))
        expected = {
            'RBP': {
                'scifi': scifi
            },
            'miRNA': {
                'fake01': fake01,
                'fake02': fake02
            },
        }

        got= utils.get_regulators()
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
                  'year': '1870'
                }
            ]
        }

        got = utils.parse_experiment(path.join(datadir, 'regulators', 'RBP', 'scifi.json'))
        self.assertEqual(expected, got)
