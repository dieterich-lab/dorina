# vim: set fileencoding=utf-8 :

import unittest
from os import path
from argparse import Namespace
from pybedtools import BedTool
from dorina import config
from dorina import utils
from dorina import run

datadir = path.join(path.dirname(path.abspath(__file__)), 'data')

class TestAnalyseWithoutOptions(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None

    def test_analyse_all_regions_seta_single(self):
        """Test run.analyse() on all regions with a single regulator"""
        expected = [
            dict(track="chr1", gene="gene01.01", data_source='scifi', score=5, site="scifi_cds",
                 location="chr1:250-260", strand="+"),
            dict(track="chr1", gene="gene01.01", data_source='scifi', score=5, site="scifi_cds",
                 location="chr1:250-260", strand="+"),
            dict(track="chr1", gene="gene01.02", data_source='scifi', score=5, site="scifi_intron",
                 location="chr1:2350-2360", strand="+")
        ]
        got = run.analyse('hg19', set_a=['scifi'], match_a='any', region_a='any', datadir=datadir)
        self.assertEqual(expected, got)

        got = run.analyse('hg19', set_a=['scifi'], datadir=datadir)
        self.assertEqual(expected, got)

    def test_analyse_CDS_seta_single(self):
        """Test run.analyse() on CDS regions with a single regulator"""
        expected = [
            dict(track="chr1", gene="gene01.01", data_source='scifi', score=5, site="scifi_cds",
                 location="chr1:250-260", strand="+")
        ]
        got = run.analyse('hg19', set_a=['scifi'], match_a='any', region_a='CDS', datadir=datadir)
        self.assertEqual(expected, got)

        got = run.analyse('hg19', set_a=['scifi'], region_a='CDS', datadir=datadir)
        self.assertEqual(expected, got)

    def test_analyse_intergenic_seta_single(self):
        """Test run.analyse() on intergenic regions with a single regulator"""
        expected = [
            dict(track="chr1", gene="intergenic01.01", data_source='scifi', score=5, site="scifi_intergenic",
                 location="chr1:1250-1260", strand=".")
        ]
        got = run.analyse('hg19', set_a=['scifi'], match_a='any', region_a='intergenic', datadir=datadir)
        self.assertEqual(expected, got)

        got = run.analyse('hg19', set_a=['scifi'], region_a='intergenic', datadir=datadir)
        self.assertEqual(expected, got)

    def test_analyse_all_regions_seta_any(self):
        """Test run.analyse() on all regions with two regulators with match to any regulator"""
        expected = [
            dict(track="chr1", gene="gene01.01", data_source='fake01', score=5, site="fake01_cds",
                 location="chr1:255-265", strand="+"),
            dict(track="chr1", gene="gene01.01", data_source='fake01', score=5, site="fake01_cds",
                 location="chr1:255-265", strand="+"),
            dict(track="chr1", gene="gene01.02", data_source='fake02', score=5, site="fake02_intron",
                 location="chr1:2450-2460", strand="+"),
            dict(track="chr1", gene="gene01.02", data_source='fake02', score=5, site="fake02_intron",
                 location="chr1:2450-2460", strand="+")
        ]
        got = run.analyse('hg19', set_a=['fake01', 'fake02'], match_a='any', region_a='any', datadir=datadir)
        self.assertEqual(expected, got)

        got = run.analyse('hg19', set_a=['fake01', 'fake02'], datadir=datadir)
        self.assertEqual(expected, got)

    def test_analyse_all_regions_seta_all(self):
        """Test run.analyse() on all regions with two regulators with match to all regulators"""
        expected = [
            dict(track="chr1", gene="gene01.01", data_source='scifi', score=5, site="scifi_cds",
                 location="chr1:250-265", strand="+"),
            dict(track="chr1", gene="gene01.01", data_source='scifi', score=5, site="scifi_cds",
                 location="chr1:250-265", strand="+")
        ]
        got = run.analyse('hg19', set_a=['scifi', 'fake01'], match_a='all', region_a='any', datadir=datadir)
        self.assertEqual(expected, got)

        got = run.analyse('hg19', set_a=['scifi', 'fake01'], match_a='all', datadir=datadir)
        self.assertEqual(expected, got)

    def test_get_genome_bedtool(self):
        """Test run._get_genome_bedtool()"""
        # should raise a ValueError for an invalid region
        self.assertRaises(ValueError, run._get_genome_bedtool, 'hg19', 'invalid', datadir)

        expected = BedTool(path.join(utils.get_genome_by_name('hg19', datadir), 'all.gff'))
        got = run._get_genome_bedtool('hg19', 'any', datadir)
        self.assertEqual(expected, got)

        expected = BedTool(path.join(utils.get_genome_by_name('hg19', datadir), 'cds.gff'))
        got = run._get_genome_bedtool('hg19', 'CDS', datadir)
        self.assertEqual(expected, got)

        expected = BedTool(path.join(utils.get_genome_by_name('hg19', datadir), '3_utr.gff'))
        got = run._get_genome_bedtool('hg19', '3prime', datadir)
        self.assertEqual(expected, got)

        expected = BedTool(path.join(utils.get_genome_by_name('hg19', datadir), '5_utr.gff'))
        got = run._get_genome_bedtool('hg19', '5prime', datadir)
        self.assertEqual(expected, got)

        expected = BedTool(path.join(utils.get_genome_by_name('hg19', datadir), 'intron.gff'))
        got = run._get_genome_bedtool('hg19', 'intron', datadir)
        self.assertEqual(expected, got)

        expected = BedTool(path.join(utils.get_genome_by_name('hg19', datadir), 'intergenic.gff'))
        got = run._get_genome_bedtool('hg19', 'intergenic', datadir)
        self.assertEqual(expected, got)


    def test_get_regulator_bedtool(self):
        """Test run._get_regulator_bedtool()"""
        expected = BedTool('%s.bed' % utils.get_regulator_by_name('scifi', datadir))
        got = run._get_regulator_bedtool('scifi', datadir)
        self.assertEqual(expected, got)


    def test_parse_results(self):
        """Test run._parse_results()"""
        results = run._analyse('hg19', set_a=['scifi', 'fake01'], match_a='all', region_a='any', datadir=datadir)
        expected = [
            dict(track="chr1", gene="gene01.01", data_source='scifi', score=5, site="scifi_cds",
                 location="chr1:250-265", strand="+"),
            dict(track="chr1", gene="gene01.01", data_source='scifi', score=5, site="scifi_cds",
                 location="chr1:250-265", strand="+")
        ]

        got = run._parse_results(results)

        self.assertEqual(expected, got)
