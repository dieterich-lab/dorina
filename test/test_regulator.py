#!/usr/bin/env python
# -*- coding: utf-8

from __future__ import unicode_literals
import unittest
import json
from os import path
from dorina import utils
from dorina.regulator import Regulator
from pybedtools import BedTool


class TestListDataWithoutOptions(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.datadir = path.join(path.dirname(path.abspath(__file__)), 'data')
        self.Regulator = Regulator.init(self.datadir)

    def TearDown(self):
        self.maxDiff = None
        self.datadir = None
        self.Regulator = None

    def test_regulator_all(self):
        """Test Regulator.all"""
        basedir = path.join(self.datadir, 'regulators', 'h_sapiens', 'hg19')

        scifi_path = path.join(self.datadir, basedir, 'PARCLIP_scifi.json')
        with open(scifi_path) as fh:
            scifi = json.load(fh)[0]
        scifi['file'] = scifi_path

        fake_path = path.join(self.datadir, basedir, 'PICTAR_fake.json')
        with open(fake_path) as fh:
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

        got = Regulator.all()
        self.maxDiff = None

        self.assertTrue("h_sapiens" in got)
        self.assertTrue("hg18" in got["h_sapiens"])
        self.assertTrue("hg19" in got["h_sapiens"])
        self.assertEqual(expected_hg19, got["h_sapiens"]["hg19"])

    def test_regulator(self):
        """Test Regulator class methods"""
        with self.assertRaises(ValueError):
            Regulator.from_name("invalid", "hg19")

        expected = path.join(self.datadir, 'regulators', 'h_sapiens', 'hg19', 'PARCLIP_scifi')
        got = Regulator.from_name("PARCLIP_scifi", "hg19").basename
        self.assertEqual(expected, got)

        expected = path.join(self.datadir, 'regulators', 'h_sapiens', 'hg19', 'PICTAR_fake')
        got = Regulator.from_name("PICTAR_fake02", "hg19").basename
        self.assertEqual(expected, got)

        expected = path.join(self.datadir, 'manual.bed')
        got = Regulator.from_name(expected).path
        self.assertEqual(expected, got)

        # Make sure that the assembly is not ignored when the
        # regulator name is not unique in the data directory.  Here we
        # have PICTAR_fake in hg18 and in hg19.
        expected = path.join(self.datadir, 'regulators', 'h_sapiens', 'hg19', 'PICTAR_fake')
        got = Regulator.from_name("PICTAR_fake01", "hg19").basename
        self.assertEqual(expected, got)

        expected = path.join(self.datadir, 'regulators', 'h_sapiens', 'hg18', 'PICTAR_fake')
        got = Regulator.from_name("PICTAR_fake01", "hg18").basename
        self.assertEqual(expected, got)

    def test_make_regulator_bed(self):
        """Test regulator.bed"""
        filename = path.join(self.datadir, 'regulators', 'h_sapiens', 'hg19', 'PARCLIP_scifi.bed')
        expected = BedTool(filename).bed6()
        got = Regulator.from_name("PARCLIP_scifi", "hg19").bed
        self.assertEqual(expected, got)

        manual = path.join(self.datadir, 'manual.bed')
        expected = BedTool(manual).bed6()
        got = Regulator.from_name(manual).bed
        self.assertEqual(expected, got)
