# vim: set fileencoding=utf-8 :

import unittest
from os import path
from argparse import Namespace
import dorina
from dorina import config

datadir = path.join(path.dirname(path.abspath(__file__)), 'data')
utils = dorina.utils.DorinaUtils(datadir)

class TestListDataWithoutOptions(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None


    def test_get_genomes(self):
        """Test utils.get_genomes()"""
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

        got = utils.get_genomes()
        self.assertEqual(expected, got)


    def test_get_regulators(self):
        """Test utils.get_regulators()"""
        basedir = path.join(datadir, 'regulators', 'h_sapiens', 'hg19')
        scifi_path = path.join(datadir, basedir, 'PARCLIP_scifi.json')
        scifi = utils.parse_experiment(scifi_path)[0]
        scifi['file'] = scifi_path
        fake_path = path.join(datadir, basedir, 'PICTAR_fake.json')
        experiments = utils.parse_experiment(fake_path)
        for exp in experiments:
            exp['file'] = fake_path

        expected_hg19 = {
            'PARCLIP_scifi': scifi,
            'PICTAR_fake01': experiments[0],
            'PICTAR_fake02': experiments[1],
            'PICTAR_fake023': experiments[2],
            'fake024|Pictar': experiments[3]
        }

        got = utils.get_regulators()
        self.maxDiff = None

        self.assertTrue("h_sapiens" in got)
        self.assertTrue("hg18" in got["h_sapiens"])
        self.assertTrue("hg19" in got["h_sapiens"])
        self.assertEqual(expected_hg19, got["h_sapiens"]["hg19"])


    def test_parse_experiment(self):
        """Test utils.parse_experiment()"""
        expected = [{
            'id': 'PARCLIP_scifi',
            'experiment': 'PARCLIP',
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
        }]

        basedir = path.join(datadir, 'regulators', 'h_sapiens', 'hg19')
        got = utils.parse_experiment(path.join(datadir, basedir, 'PARCLIP_scifi.json'))
        self.assertEqual(expected, got)


    def test_get_genome_by_name(self):
        """Test utils.get_genome_by_name()"""
        got = utils.get_genome_by_name("invalid")
        self.assertIsNone(got)

        expected = path.join(datadir, 'genomes', 'h_sapiens', 'hg19')
        got = utils.get_genome_by_name("hg19")
        self.assertEqual(expected, got)


    def test_get_regulator_by_name(self):
        """Test utils.get_regulator_by_name()"""
        got = utils.get_regulator_by_name("invalid")
        self.assertIsNone(got)

        expected = path.join(datadir, 'regulators', 'h_sapiens', 'hg19', 'PARCLIP_scifi')
        got = utils.get_regulator_by_name("PARCLIP_scifi")
        self.assertEqual(expected, got)

        expected = path.join(datadir, 'regulators', 'h_sapiens', 'hg19', 'PICTAR_fake')
        got = utils.get_regulator_by_name("PICTAR_fake02")
        self.assertEqual(expected, got)

        expected = path.join(datadir, 'manual')
        got = utils.get_regulator_by_name(expected)
        self.assertEqual(expected, got)

        # Make sure that the assembly is not ignored when the
        # regulator name is not unique in the data directory.  Here we
        # have PICTAR_fake in hg18 and in hg19.
        expected = path.join(datadir, 'regulators', 'h_sapiens', 'hg19', 'PICTAR_fake')
        got = utils.get_regulator_by_name("PICTAR_fake")
        self.assertEqual(expected, got)

    def test_get_genes(self):
        """Test utils.get_genes()"""
        expected = ['gene01.01', 'gene01.02']
        got = utils.get_genes('hg19')
        self.assertEqual(expected, got)

        expected = []
        got = utils.get_genes('invalid')
        self.assertEqual(expected, got)
