import unittest

from pandas import DataFrame
from pandas.testing import assert_frame_equal

from dorina.expression import (retrieve_study, retrieve_fpkm_from_study,
                               calculate_expressed_genes)

try:
    from unittest import mock
except ImportError:
    import mock


class TestExpression(unittest.TestCase):
    def setUp(self):
        self.retrieve_study = retrieve_study
        self.fpkm_from_study = retrieve_fpkm_from_study
        self.calculate = calculate_expressed_genes
        self.get_run_by_organism_response = [
            {'ASSEMBLY_USED': 'test',
             'BEDGRAPH_LOCATION': 'ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/DRR000/DRR000897/DRR000897.bedgraph',
             'BIGWIG_LOCATION': 'ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/DRR000/DRR000897/DRR000897.bw',
             'BIOREP_ID': 'DRR000897',
             'CRAM_LOCATION': 'ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/DRR000/DRR000897/DRR000897.cram',
             'ENA_LAST_UPDATED': 'Fri Jun 19 2015 17:39:46',
             'LAST_PROCESSED_DATE': 'Sun Sep 20 2015 01:38:35',
             'MAPPING_QUALITY': 98,
             'ORGANISM': 'homo_sapiens',
             'REFERENCE_ORGANISM': 'homo_sapiens',
             'RUN_IDS': 'DRR000897',
             'SAMPLE_IDS': 'SAMD00004140',
             'STATUS': 'Complete',
             'STUDY_ID': 'DRP000366'},
            {'ASSEMBLY_USED': 'GRCh38',
             'BEDGRAPH_LOCATION': 'ftp://DRR001176.bedgraph',
             'BIGWIG_LOCATION': 'ftp://DRR001176.bw',
             'BIOREP_ID': 'DRR001176',
             'CRAM_LOCATION': 'ftp://.cram',
             'ENA_LAST_UPDATED': 'Fri Jun 19 2015 17:39:49',
             'LAST_PROCESSED_DATE': 'Sun Aug 30 2015 03:58:31',
             'MAPPING_QUALITY': 97,
             'ORGANISM': 'homo_sapiens',
             'REFERENCE_ORGANISM': 'homo_sapiens',
             'RUN_IDS': 'DRR001176',
             'SAMPLE_IDS': 'SAMD00004026',
             'STATUS': 'Fastq_failed',
             'STUDY_ID': 'DRP000425'},
            {'ASSEMBLY_USED': 'GRCh38',
             'BEDGRAPH_LOCATION': 'ftp://DRR001174.bedgraph',
             'BIGWIG_LOCATION': 'ftp://DRR001174.bw',
             'BIOREP_ID': 'DRR001174',
             'CRAM_LOCATION': 'ftp://DRR001174.cram',
             'ENA_LAST_UPDATED': 'Fri Jun 19 2015 17:39:49',
             'LAST_PROCESSED_DATE': 'Fri Aug 28 2015 04:36:31',
             'MAPPING_QUALITY': 98,
             'ORGANISM': 'homo_sapiens',
             'REFERENCE_ORGANISM': 'homo_sapiens',
             'RUN_IDS': 'DRR001174',
             'SAMPLE_IDS': 'SAMD00004025',
             'STATUS': 'Complete',
             'STUDY_ID': 'DRP003703'}]
        self.get_studies_by_organism_response = [
            {'ASSEMBLY_USED': 'GRCh38',
             'EXONS_FPKM_COUNTS_FTP_LOCATION': 'ftp://ena/DRP003703/homo_sapiens/exons.fpkm.dexseq.tsv',
             'EXONS_RAW_COUNTS_FTP_LOCATION': 'ftp://ena/DRP003703/homo_sapiens/exons.raw.dexseq.tsv',
             'EXONS_TPM_COUNTS_FTP_LOCATION': 'ftp://ena/DRP003703/homo_sapiens/exons.tpm.dexseq.tsv',
             'GENES_FPKM_COUNTS_FTP_LOCATION': 'ftp://ena/DRP003703/homo_sapiens/genes.fpkm.htseq2.tsv',
             'GENES_RAW_COUNTS_FTP_LOCATION': 'ftp://ena/DRP003703/homo_sapiens/genes.raw.htseq2.tsv',
             'GENES_TPM_COUNTS_FTP_LOCATION': 'ftp://ena/DRP003703/homo_sapiens/genes.tpm.htseq2.tsv',
             'GTF_USED': 'Homo_sapiens.GRCh38.79.gtf.gz',
             'LAST_PROCESSED_DATE': 'Fri Jan 12 2018 04:13:45',
             'ORGANISM': 'homo_sapiens',
             'REFERENCE_ORGANISM': 'homo_sapiens',
             'SOFTWARE_VERSIONS_FTP_LOCATION': 'ftp://ena/DRP003703/homo_sapiens/irap.versions.tsv',
             'STATUS': 'Complete',
             'STUDY_ID': 'DRP003703'}]
        self.study_gene_fpkm_response = DataFrame(
            {'ENSG00000000003': {'DRR078748': 0.0,
                                 'DRR078749': 0.0,
                                 'DRR078750': 0.0,
                                 'DRR078751': 3.8,
                                 'DRR078752': 2.23,
                                 'DRR078753': 0.86,
                                 'DRR078754': 0.0,
                                 'DRR078755': 0.0,
                                 'DRR078757': 0.33,
                                 'DRR078758': 0.05,
                                 'DRR078759': 0.38,
                                 'DRR078760': 12.65},
             'ENSG00000000005': {'DRR078748': 0.0,
                                 'DRR078749': 0.0,
                                 'DRR078750': 0.0,
                                 'DRR078751': 0.0,
                                 'DRR078752': 0.0,
                                 'DRR078753': 0.0,
                                 'DRR078754': 0.0,
                                 'DRR078755': 0.0,
                                 'DRR078757': 0.0,
                                 'DRR078758': 0.0,
                                 'DRR078759': 0.0,
                                 'DRR078760': 0.0},
             'ENSG00000000419': {'DRR078748': 0.0,
                                 'DRR078749': 0.0,
                                 'DRR078750': 0.0,
                                 'DRR078751': 0.0,
                                 'DRR078752': 0.0,
                                 'DRR078753': 0.0,
                                 'DRR078754': 0.0,
                                 'DRR078755': 0.0,
                                 'DRR078757': 0.0,
                                 'DRR078758': 0.0,
                                 'DRR078759': 0.0,
                                 'DRR078760': 0.0},
             'ENSG00000000457': {'DRR078748': 0.31,
                                 'DRR078749': 0.48,
                                 'DRR078750': 0.0,
                                 'DRR078751': 0.83,
                                 'DRR078752': 0.0,
                                 'DRR078753': 0.0,
                                 'DRR078754': 0.0,
                                 'DRR078755': 0.3,
                                 'DRR078757': 0.22,
                                 'DRR078758': 0.21,
                                 'DRR078759': 0.0,
                                 'DRR078760': 0.0},
             'ENSG00000000460': {'DRR078748': 0.0,
                                 'DRR078749': 0.28,
                                 'DRR078750': 0.0,
                                 'DRR078751': 0.61,
                                 'DRR078752': 0.14,
                                 'DRR078753': 0.0,
                                 'DRR078754': 0.0,
                                 'DRR078755': 0.0,
                                 'DRR078757': 0.75,
                                 'DRR078758': 0.48,
                                 'DRR078759': 0.97,
                                 'DRR078760': 0.23}})

    def TearDown(self):
        self.retrieve_study = None
        self.get_run_by_organism_response = None
        self.get_studies_by_organism_response = None
        self.study_gene_fpkm = None
        self.study_gene_fpkm_response = None

    def test_bad_assembly(self):
        with self.assertRaises(KeyError):
            self.retrieve_study('test')

    @mock.patch('dorina.expression.retrieve')
    def test_bad_condition(self, mock_retrieve):
        mock_retrieve.get_run_by_organism.return_value = []
        with self.assertRaises(ValueError):
            self.retrieve_study('hg19', 'test')

    @mock.patch('dorina.expression.retrieve')
    def test_retrieve_study(self, mock_retrieve):
        mock_retrieve.get_run_by_organism.return_value = \
            self.get_studies_by_organism_response

        mock_retrieve.get_studies_by_organism.return_value = \
            self.get_run_by_organism_response

        self.assertEqual(self.retrieve_study('GRCh38'),
                         [self.get_run_by_organism_response[2]])
        mock_retrieve.get_run_by_organism.assert_called_once_with(
            organism='homo_sapiens', condition=None)

        mock_retrieve.get_studies_by_organism.assert_called_once_with(
            'homo_sapiens')

    def test_retrieve_fpkm_from_study_no_studies(self):
        with self.assertRaises(ValueError):
            self.fpkm_from_study([])

    @mock.patch('dorina.expression.read_table')
    def test_retrieve_fpkm_from_study(self, mock_read_table):
        mock_read_table.return_value = self.study_gene_fpkm_response
        data = self.fpkm_from_study(
            [{'GENES_FPKM_COUNTS_FTP_LOCATION': 'test'}])
        assert_frame_equal(data, self.study_gene_fpkm_response)

    @mock.patch('dorina.expression.read_table')
    def test_calculate_expressed_genes(self, mock_read_table):
        mock_read_table.return_value = self.study_gene_fpkm_response.T
        data = self.fpkm_from_study(
            [{'GENES_FPKM_COUNTS_FTP_LOCATION': 'test'}])
        self.assertEqual(self.calculate(data), ['ENSG00000000003'])
        self.assertEqual(
            sorted(self.calculate(data, -1)),
            sorted(['ENSG00000000003', 'ENSG00000000005', 'ENSG00000000419',
                    'ENSG00000000457', 'ENSG00000000460']))
        self.assertEqual(self.calculate(data, 200), [])


if __name__ == '__main__':
    import nose

    nose.runmodule()
