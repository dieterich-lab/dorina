# vim: set fileencoding=utf-8 :

import unittest
from os import path
from dorina import utils
from dorina.regulator import Regulator
from pybedtools import BedTool

datadir = path.join(path.dirname(path.abspath(__file__)), 'data')
utils = utils.DorinaUtils(datadir)

class TestListDataWithoutOptions(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None

    def test_regulator(self):
        """Test Regulator class methods"""
        # TODO: expect an error, not None!
        #got = utils.make_regulator("invalid", "hg19")
        #self.assertIsNone(got)

        expected = path.join(datadir, 'regulators', 'h_sapiens', 'hg19', 'PARCLIP_scifi')
        got = Regulator.from_name(utils.regulators, "PARCLIP_scifi", "hg19").basename
        self.assertEqual(expected, got)

        expected = path.join(datadir, 'regulators', 'h_sapiens', 'hg19', 'PICTAR_fake')
        got = Regulator.from_name(utils.regulators, "PICTAR_fake02", "hg19").basename
        self.assertEqual(expected, got)

        expected = path.join(datadir, 'manual.bed')
        got = Regulator.from_name(utils.regulators, expected).path
        self.assertEqual(expected, got)

        # Make sure that the assembly is not ignored when the
        # regulator name is not unique in the data directory.  Here we
        # have PICTAR_fake in hg18 and in hg19.
        expected = path.join(datadir, 'regulators', 'h_sapiens', 'hg19', 'PICTAR_fake')
        got = Regulator.from_name(utils.regulators, "PICTAR_fake01", "hg19").basename
        self.assertEqual(expected, got)

        expected = path.join(datadir, 'regulators', 'h_sapiens', 'hg18', 'PICTAR_fake')
        got = Regulator.from_name(utils.regulators, "PICTAR_fake01", "hg18").basename
        self.assertEqual(expected, got)


    def test_make_regulator_bed(self):
        """Test regulator.bed"""
        filename = path.join(datadir, 'regulators', 'h_sapiens', 'hg19', 'PARCLIP_scifi.bed')
        expected = BedTool(filename).bed6()
        got = Regulator.from_name(utils.regulators, "PARCLIP_scifi", "hg19").bed
        self.assertEqual(expected, got)

        manual = path.join(datadir, 'manual.bed')
        expected = BedTool(manual).bed6()
        got = Regulator.from_name(utils.regulators, manual).bed
        self.assertEqual(expected, got)
