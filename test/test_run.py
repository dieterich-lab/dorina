# vim: set fileencoding=utf-8 :

import unittest
from os import path
from argparse import Namespace
from pybedtools import BedTool

from dorina import config
from dorina import run

datadir = path.join(path.dirname(path.abspath(__file__)), 'data')
run = run.Dorina(datadir)
utils = run.utils

class TestAnalyseWithoutOptions(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None

    def test_analyse_all_regions_seta_single(self):
        """Test run.analyse() on all regions with a single regulator"""
        bed_str = """chr1   doRiNA2 gene    1   1000    .   +   .   ID=gene01.01    chr1    250 260 PARCLIP#scifi*scifi_cds 5   +
        chr1    doRiNA2 gene    2001    3000    .   +   .   ID=gene01.02    chr1    2350    2360    PARCLIP#scifi*scifi_intron  5   +"""
        expected = BedTool(bed_str, from_string=True)
        got = run.analyse('hg19', set_a=['PARCLIP_scifi'], match_a='any', region_a='any')
        self.assertEqual(expected, got)

        got = run.analyse('hg19', set_a=['PARCLIP_scifi'])
        self.assertEqual(expected, got)

    def test_analyse_CDS_seta_single(self):
        """Test run.analyse() on CDS regions with a single regulator"""
        bed_str = """chr1   doRiNA2 CDS 201 300 .   +   0   ID=gene01.01    chr1    250 260 PARCLIP#scifi*scifi_cds 5   +"""
        expected = BedTool(bed_str, from_string=True)
        got = run.analyse('hg19', set_a=['PARCLIP_scifi'], match_a='any', region_a='CDS')
        self.assertEqual(expected, got)

        got = run.analyse('hg19', set_a=['PARCLIP_scifi'], region_a='CDS')
        self.assertEqual(expected, got)

    def test_analyse_intergenic_seta_single(self):
        """Test run.analyse() on intergenic regions with a single regulator"""
        bed_str = """chr1   doRiNA2 intergenic  1001    2000    .   .   .   ID=intergenic01.01  chr1    1250    1260    PARCLIP#scifi*scifi_intergenic  5   ."""
        expected = BedTool(bed_str, from_string=True)
        got = run.analyse('hg19', set_a=['PARCLIP_scifi'], match_a='any', region_a='intergenic')
        self.assertEqual(expected, got)

        got = run.analyse('hg19', set_a=['PARCLIP_scifi'], region_a='intergenic')
        self.assertEqual(expected, got)

    def test_analyse_all_regions_seta_any(self):
        """Test run.analyse() on all regions with two regulators with match to any regulator"""
        bed_str = """chr1   doRiNA2 gene    1   1000    .   +   .   ID=gene01.01    chr1    255 265 PICTAR#fake01*fake01_cds    5   +
        chr1    doRiNA2 gene    2001    3000    .   +   .   ID=gene01.02    chr1    2450    2460    PICTAR#fake02*fake02_intron 5   +"""
        expected = BedTool(bed_str, from_string=True)
        got = run.analyse('hg19', set_a=['PICTAR_fake01', 'PICTAR_fake02'], match_a='any', region_a='any')
        self.assertMultiLineEqual(str(expected), str(got))

        got = run.analyse('hg19', set_a=['PICTAR_fake01', 'PICTAR_fake02'])
        self.assertEqual(expected, got)

    def test_analyse_all_regions_seta_all(self):
        """Test run.analyse() on all regions with two regulators with match to all regulators"""
        bed_str = """chr1   doRiNA2 gene    1   1000    .   +   .   ID=gene01.01    chr1    250 260 PARCLIP#scifi*scifi_cds 5   +
        chr1    doRiNA2 gene    1   1000    .   +   .   ID=gene01.01    chr1    255 265 PICTAR#fake01*fake01_cds    5   +"""
        expected = BedTool(bed_str, from_string=True)
        got = run.analyse('hg19', set_a=['PARCLIP_scifi', 'PICTAR_fake01'], match_a='all', region_a='any')
        self.assertMultiLineEqual(str(expected), str(got))

        got = run.analyse('hg19', set_a=['PARCLIP_scifi', 'PICTAR_fake01'], match_a='all')
        self.assertEqual(expected, got)

    def test_analyse_all_regions_seta_and_setb(self):
        """Test run.analyse() on all regions with any regulator from set A and any regulator from set B matching"""
        bed_str = """chr1   doRiNA2 gene    1   1000    .   +   .   ID=gene01.01    chr1    250 260 PARCLIP#scifi*scifi_cds 5   +
        chr1    doRiNA2 gene    1   1000    .   +   .   ID=gene01.01    chr1    255 265 PICTAR#fake01*fake01_cds    5   +"""
        expected = BedTool(bed_str, from_string=True)
        got = run.analyse('hg19', set_a=['PARCLIP_scifi'], match_a='any', region_a='any',
                          set_b=['PICTAR_fake01'], match_b='any', region_b='any',
                          combine='and')
        self.assertMultiLineEqual(str(expected), str(got))

    def test_analyse_all_regions_seta_or_setb(self):
        """Test run.analyse() on all regions with any regulator from set A or any regulator from set B matching"""
        bed_str = """chr1   doRiNA2 gene    1   1000    .   +   .   ID=gene01.01    chr1    250 260 PARCLIP#scifi*scifi_cds 5   +
        chr1    doRiNA2 gene    1   1000    .   +   .   ID=gene01.01    chr1    255 265 PICTAR#fake01*fake01_cds    5   +
        chr1    doRiNA2 gene    2001    3000    .   +   .   ID=gene01.02    chr1    2350    2360    PARCLIP#scifi*scifi_intron  5   +
        chr1    doRiNA2 gene    1   1000    .   +   .   ID=gene01.01    chr1    250 260 PARCLIP#scifi*scifi_cds 5   +
        chr1    doRiNA2 gene    1   1000    .   +   .   ID=gene01.01    chr1    255 265 PICTAR#fake01*fake01_cds    5   +"""
        expected = BedTool(bed_str, from_string=True)
        got = run.analyse('hg19', set_a=['PARCLIP_scifi'], match_a='any', region_a='any',
                          set_b=['PICTAR_fake01'], match_b='any', region_b='any',
                          combine='or')
        self.assertMultiLineEqual(str(expected), str(got))

    def test_analyse_all_regions_seta_xor_setb(self):
        """Test run.analyse() on all regions with any regulator from set A XOR any regulator from set B matching"""
        bed_str = """chr1   doRiNA2 gene    2001    3000    .   +   .   ID=gene01.02    chr1    2350    2360    PARCLIP#scifi*scifi_intron  5   +"""
        expected = BedTool(bed_str, from_string=True)
        got = run.analyse('hg19', set_a=['PARCLIP_scifi'], match_a='any', region_a='any',
                          set_b=['PICTAR_fake01'], match_b='any', region_b='any',
                          combine='xor')
        self.assertMultiLineEqual(str(expected), str(got))

    def test_analyse_all_regions_seta_not_setb(self):
        """Test run.analyse() on all regions with any regulator from set A but no regulator from set B matching"""
        bed_str = """chr1   doRiNA2 gene    2001    3000    .   +   .   ID=gene01.02    chr1    2350    2360    PARCLIP#scifi*scifi_intron  5   +"""
        expected = BedTool(bed_str, from_string=True)
        got = run.analyse('hg19', set_a=['PARCLIP_scifi'], match_a='any', region_a='any',
                          set_b=['PICTAR_fake01'], match_b='any', region_b='any',
                          combine='not')
        self.assertMultiLineEqual(str(expected), str(got))


    def test_analyse_all_regions_seta_windowed(self):
        """Test run.analyse() on all regions with all regulators from set A matching in an overlapping window"""
        bed_str = """chr1	doRiNA2	gene	250	260 .	+	.	ID=gene01.01	chr1	250	260	PARCLIP#scifi*scifi_cds	5	+
chr1	doRiNA2	gene	250	260	.	+	.	ID=gene01.01	chr1	255	265	PICTAR#fake01*fake01_cds	5	+"""
        expected = BedTool(bed_str, from_string=True)
        got = run.analyse('hg19', set_a=['PARCLIP_scifi', 'PICTAR_fake01'], match_a='all', region_a='any',
                          window_a=0)
        self.assertMultiLineEqual(str(expected), str(got))


    def test_analyse_all_regions_seta_windowed_slop(self):
        """Test run.analyse() on all regions with all regulators from set A matching in an overlapping window with slop"""
        bed_str = """chr1	doRiNA2	gene    1	1260	.	+	.	ID=gene01.01	chr1	250	260	PARCLIP#scifi*scifi_cds	5	+
chr1	doRiNA2	gene	1	1260	.	+	.	ID=gene01.01	chr1	1250	1260	PARCLIP#scifi*scifi_intergenic	5	.
chr1	doRiNA2	gene	1	1260	.	+	.	ID=gene01.01	chr1	255	265	PICTAR#fake01*fake01_cds	5	+
chr1	doRiNA2	gene	1350	3360	.	+	.	ID=gene01.02	chr1	2350	2360	PARCLIP#scifi*scifi_intron	5	+
chr1	doRiNA2	gene	1350	3360	.	+	.	ID=gene01.02	chr1	1350	1360	PICTAR#fake01*fake01_intergenic	5	."""
        expected = BedTool(bed_str, from_string=True)
        got = run.analyse('hg19', set_a=['PARCLIP_scifi', 'PICTAR_fake01'], match_a='all', region_a='any',
                          window_a=1000)
        self.assertMultiLineEqual(str(expected), str(got))


    def test_add_slop(self):
        """Test run._add_slop()"""
        slop_string = """chr1   0   560 PARCLIP#scifi*scifi_cds 5   +
        chr1    950 1560    PARCLIP#scifi*scifi_intergenic  5   .
        chr1    2050    2660    PARCLIP#scifi*scifi_intron  5   +
"""
        expected = BedTool(slop_string, from_string=True)
        got = run._add_slop(BedTool(run._get_regulator_bedtool('PARCLIP_scifi')),
                            'hg19', 300)
        self.assertEqual(expected, got)


    def test_get_genome_bedtool(self):
        """Test run._get_genome_bedtool()"""
        # should raise a ValueError for an invalid region
        self.assertRaises(ValueError, run._get_genome_bedtool, 'hg19', 'invalid')

        expected = BedTool(path.join(utils.get_genome_by_name('hg19'), 'all.gff'))
        got = run._get_genome_bedtool('hg19', 'any')
        self.assertEqual(expected, got)

        expected = BedTool(path.join(utils.get_genome_by_name('hg19'), 'cds.gff'))
        got = run._get_genome_bedtool('hg19', 'CDS')
        self.assertEqual(expected, got)

        expected = BedTool(path.join(utils.get_genome_by_name('hg19'), '3_utr.gff'))
        got = run._get_genome_bedtool('hg19', '3prime')
        self.assertEqual(expected, got)

        expected = BedTool(path.join(utils.get_genome_by_name('hg19'), '5_utr.gff'))
        got = run._get_genome_bedtool('hg19', '5prime')
        self.assertEqual(expected, got)

        expected = BedTool(path.join(utils.get_genome_by_name('hg19'), 'intron.gff'))
        got = run._get_genome_bedtool('hg19', 'intron')
        self.assertEqual(expected, got)

        expected = BedTool(path.join(utils.get_genome_by_name('hg19'), 'intergenic.gff'))
        got = run._get_genome_bedtool('hg19', 'intergenic')
        self.assertEqual(expected, got)

        expected = BedTool(path.join(utils.get_genome_by_name('hg19'), 'all.gff')).filter(
                lambda x: x.name == "gene01.02").saveas()
        got = run._get_genome_bedtool('hg19', 'any', genes=['gene01.02'])
        self.assertEqual(expected, got)


    def test_get_regulator_bedtool(self):
        """Test run._get_regulator_bedtool()"""
        expected = BedTool('%s.bed' % utils.get_regulator_by_name('PARCLIP_scifi')).bed6()
        got = run._get_regulator_bedtool('PARCLIP_scifi')
        self.assertEqual(expected, got)

        manual = path.join(datadir, 'manual.bed')
        expected = BedTool(manual).bed6()
        got = run._get_regulator_bedtool(manual)
        self.assertEqual(expected, got)
