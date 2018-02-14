import unittest

from dorina.expression import retrieve_expressed_genes

try:
    from unittest import mock
except ImportError:
    import mock

from nose.tools import assert_raises


class TestExpression(unittest.TestCase):
    def setUp(self):
        self.retrieve_expressed_genes = retrieve_expressed_genes
        self.retrieve_get_run_by_organism_response = [
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
             'BEDGRAPH_LOCATION': 'ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/DRR001/DRR001176/DRR001176.bedgraph',
             'BIGWIG_LOCATION': 'ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/DRR001/DRR001176/DRR001176.bw',
             'BIOREP_ID': 'DRR001176',
             'CRAM_LOCATION': 'ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/DRR001/DRR001176/DRR001176.cram',
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
             'BEDGRAPH_LOCATION': 'ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/DRR001/DRR001174/DRR001174.bedgraph',
             'BIGWIG_LOCATION': 'ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/DRR001/DRR001174/DRR001174.bw',
             'BIOREP_ID': 'DRR001174',
             'CRAM_LOCATION': 'ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/DRR001/DRR001174/DRR001174.cram',
             'ENA_LAST_UPDATED': 'Fri Jun 19 2015 17:39:49',
             'LAST_PROCESSED_DATE': 'Fri Aug 28 2015 04:36:31',
             'MAPPING_QUALITY': 98,
             'ORGANISM': 'homo_sapiens',
             'REFERENCE_ORGANISM': 'homo_sapiens',
             'RUN_IDS': 'DRR001174',
             'SAMPLE_IDS': 'SAMD00004025',
             'STATUS': 'Complete',
             'STUDY_ID': 'DRP003703'}]
        self.retrieve_get_studies_by_organism = [{'ASSEMBLY_USED': 'GRCh38',
                                                  'EXONS_FPKM_COUNTS_FTP_LOCATION': 'ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/studies/ena/DRP003703/homo_sapiens/exons.fpkm.dexseq.tsv',
                                                  'EXONS_RAW_COUNTS_FTP_LOCATION': 'ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/studies/ena/DRP003703/homo_sapiens/exons.raw.dexseq.tsv',
                                                  'EXONS_TPM_COUNTS_FTP_LOCATION': 'ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/studies/ena/DRP003703/homo_sapiens/exons.tpm.dexseq.tsv',
                                                  'GENES_FPKM_COUNTS_FTP_LOCATION': 'ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/studies/ena/DRP003703/homo_sapiens/genes.fpkm.htseq2.tsv',
                                                  'GENES_RAW_COUNTS_FTP_LOCATION': 'ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/studies/ena/DRP003703/homo_sapiens/genes.raw.htseq2.tsv',
                                                  'GENES_TPM_COUNTS_FTP_LOCATION': 'ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/studies/ena/DRP003703/homo_sapiens/genes.tpm.htseq2.tsv',
                                                  'GTF_USED': 'Homo_sapiens.GRCh38.79.gtf.gz',
                                                  'LAST_PROCESSED_DATE': 'Fri Jan 12 2018 04:13:45',
                                                  'ORGANISM': 'homo_sapiens',
                                                  'REFERENCE_ORGANISM': 'homo_sapiens',
                                                  'SOFTWARE_VERSIONS_FTP_LOCATION': 'ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/studies/ena/DRP003703/homo_sapiens/irap.versions.tsv',
                                                  'STATUS': 'Complete',
                                                  'STUDY_ID': 'DRP003703'}]
        self.study_gene_fpkm = {
            'DRR078748': {0: 0., 1: 0.0, 2: 0.0, 3: 0.31, 4: 0.0},
            'DRR078749': {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.48, 4: 0.28},
            'DRR078750': {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0},
            'DRR078751': {0: 3.8, 1: 0.0, 2: 0.0, 3: 0.83, 4: 0.61},
            'DRR078752': {0: 2.23, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.14},
            'DRR078753': {0: 0.86, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0},
            'DRR078754': {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0},
            'DRR078755': {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.3, 4: 0.0},
            'DRR078757': {0: 0.33, 1: 0.0, 2: 0.0, 3: 0.22, 4: 0.75},
            'DRR078758': {0: 0.05, 1: 0.0, 2: 0.0, 3: 0.21, 4: 0.48},
            'DRR078759': {0: 0.38, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.97},
            'DRR078760': {0: 12.65, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.23},
            'Gene ID': {0: 'ENSG00000000003',
                        1: 'ENSG00000000005',
                        2: 'ENSG00000000419',
                        3: 'ENSG00000000457',
                        4: 'ENSG00000000460'}}

    def TearDown(self):
        self.retrieve_expressed_genes = None
        self.retrieve_get_run_by_organism_response = None
        self.retrieve_get_studies_by_organism = None
        self.study_gene_fpkm = None

    def test_bad_assembly(self):
        assert_raises(KeyError, self.retrieve_expressed_genes, 'test')

    @mock.patch('dorina.expression.retrieve')
    def test_no_condition_found(self, mock_retrieve):
        retrieve = mock_retrieve()
        retrieve.get_run_by_organism.return_value = []
        assert_raises(ValueError, self.retrieve_expressed_genes, 'hg19', 'test')

    @mock.patch('dorina.expression.retrieve')
    def test_no_studies_after_filter(self, mock_retrieve):
        retrieve = mock_retrieve()
        retrieve.get_run_by_organism.return_value = \
            list(self.retrieve_get_run_by_organism_response[0])
        assert_raises(ValueError, self.retrieve_expressed_genes, 'hg19')

    # @mock.patch('dorina.expression.retrieve')
    # @mock.patch('dorina.expression.read_table')
    # def test_retrieve_for_organism(self, mock_retrieve, mock_read_table):
    #     mock_retrieve.get_run_by_organism().return_value = \
    #         self.retrieve_get_run_by_organism_response
    #     mock_retrieve.get_studies_by_organism().return_value = \
    #         self.retrieve_get_studies_by_organism
    #     mock_read_table.return_value = pd.DataFrame(self.study_gene_fpkm)
    #     assert self.retrieve_expressed_genes('hg19') == ['ENSG00000000003']


if __name__ == '__main__':
    import nose

    nose.runmodule()
