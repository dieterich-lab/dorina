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
        expected = """chr1	doRiNA2	gene	1	1000	.	+	.	gene01.01	chr1	250	260	scifi_cds	5	+	250	260
chr1	doRiNA2	CDS	201	300	.	+	0	gene01.01	chr1	250	260	scifi_cds	5	+	250	260
chr1	doRiNA2	gene	2001	3000	.	+	.	gene01.02	chr1	2350	2360	scifi_intron	5	+	2350	2360
"""
        got = run.analyse('hg19', set_a=['scifi'], match_a='any', region_a='any', datadir=datadir)
        self.assertMultiLineEqual(expected, str(got))

        got = run.analyse('hg19', set_a=['scifi'], datadir=datadir)
        self.assertMultiLineEqual(expected, str(got))

    def test_analyse_CDS_seta_single(self):
        """Test run.analyse() on CDS regions with a single regulator"""
        expected = """chr1	doRiNA2	CDS	201	300	.	+	0	gene01.01	chr1	250	260	scifi_cds	5	+	250	260
"""
        got = run.analyse('hg19', set_a=['scifi'], match_a='any', region_a='CDS', datadir=datadir)
        self.assertMultiLineEqual(expected, str(got))

        got = run.analyse('hg19', set_a=['scifi'], region_a='CDS', datadir=datadir)
        self.assertMultiLineEqual(expected, str(got))

    def test_analyse_intergenic_seta_single(self):
        """Test run.analyse() on intergenic regions with a single regulator"""
        expected = """chr1	doRiNA2	intergenic	1001	2000	.	.	.	intergenic01.01	chr1	1250	1260	scifi_intergenic	5		1250	1260
"""
        got = run.analyse('hg19', set_a=['scifi'], match_a='any', region_a='intergenic', datadir=datadir)
        self.assertMultiLineEqual(expected, str(got))

        got = run.analyse('hg19', set_a=['scifi'], region_a='intergenic', datadir=datadir)
        self.assertMultiLineEqual(expected, str(got))

    def test_analyse_all_regions_seta_any(self):
        """Test run.analyse() on all regions with two regulators with match to any regulator"""
        expected = """chr1	doRiNA2	gene	1	1000	.	+	.	gene01.01	chr1	255	265	fake01_cds	5	+	255	265
chr1	doRiNA2	CDS	201	300	.	+	0	gene01.01	chr1	255	265	fake01_cds	5	+	255	265
chr1	doRiNA2	gene	2001	3000	.	+	.	gene01.02	chr1	2450	2460	fake02_intron	5	+	2450	2460
chr1	doRiNA2	CDS	2401	2700	.	+	1	gene01.02	chr1	2450	2460	fake02_intron	5	+	2450	2460
"""
        got = run.analyse('hg19', set_a=['fake01', 'fake02'], match_a='any', region_a='any', datadir=datadir)
        self.assertMultiLineEqual(expected, str(got))

        got = run.analyse('hg19', set_a=['fake01', 'fake02'], datadir=datadir)
        self.assertMultiLineEqual(expected, str(got))

    def test_analyse_all_regions_seta_all(self):
        """Test run.analyse() on all regions with two regulators with match to all regulators"""
        expected = """chr1	doRiNA2	gene	1	1000	.	+	.	gene01.01	chr1	250	260	scifi_cds	5	+	250	260	chr1	255	265	fake01_cds	5	+	255	265
chr1	doRiNA2	CDS	201	300	.	+	0	gene01.01	chr1	250	260	scifi_cds	5	+	250	260	chr1	255	265	fake01_cds	5	+	255	265
"""
        got = run.analyse('hg19', set_a=['scifi', 'fake01'], match_a='all', region_a='any', datadir=datadir)
        self.assertMultiLineEqual(expected, str(got))

        got = run.analyse('hg19', set_a=['scifi', 'fake01'], match_a='all', datadir=datadir)
        self.assertMultiLineEqual(expected, str(got))

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
