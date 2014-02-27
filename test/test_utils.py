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
