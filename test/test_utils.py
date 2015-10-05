# vim: set fileencoding=utf-8 :

import unittest
from os import path
from argparse import Namespace
import json

import dorina
from dorina import config
from dorina.regulator import Regulator
from pybedtools import BedTool

datadir = path.join(path.dirname(path.abspath(__file__)), 'data')
utils = dorina.utils.DorinaUtils(datadir)

class TestListDataWithoutOptions(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None


    def test_genomes(self):
        """Test utils.genomes"""
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

        got = utils.genomes
        self.assertEqual(expected, got)


    def test_regulators(self):
        """Test utils.regulators"""
        basedir = path.join(datadir, 'regulators', 'h_sapiens', 'hg19')

        scifi_path = path.join(datadir, basedir, 'PARCLIP_scifi.json')
        with open(scifi_path, 'r') as fh:
            scifi = json.load(fh)[0]
        scifi['file'] = scifi_path

        fake_path = path.join(datadir, basedir, 'PICTAR_fake.json')
        with open(fake_path, 'r') as fh:
            experiments = json.load(fh)

        for exp in experiments:
            exp['file'] = fake_path

        expected_hg19 = {
            'PARCLIP_scifi': scifi,
            'PICTAR_fake01': experiments[0],
            'PICTAR_fake02': experiments[1],
            'PICTAR_fake023': experiments[2],
            'fake024|Pictar': experiments[3]
        }

        got = utils.regulators
        self.maxDiff = None

        self.assertTrue("h_sapiens" in got)
        self.assertTrue("hg18" in got["h_sapiens"])
        self.assertTrue("hg19" in got["h_sapiens"])
        self.assertEqual(expected_hg19, got["h_sapiens"]["hg19"])
