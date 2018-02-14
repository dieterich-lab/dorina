#!/usr/bin/env python
# -*- coding: utf-8
"""
Created on 09:21 14/02/2018 2018 

"""
import logging

from pandas import read_table, concat
from bioservices import RNASEQ_EBI

log = logging.getLogger(__name__)

assembly_to_organism = {'ce6': 'caenorhabditis_elegans',
                        'dm3': "drosophila_melanogaster",
                        'hg19': 'homo_sapiens', 'mm9': 'mus_musculus',
                        'GRCh38': 'homo_sapiens'}

retrieve = RNASEQ_EBI()


def retrieve_expressed_genes(assembly, condition=None, fpkm_cutoff=1):
    """
    Retrieves genes expressed expressed (fpkm > `fpkm_cutoff`) for a given
     `assembly`. Optionally filter for a experimental condition, such as 'heart'

    :param str assembly: Genomic assembly studied
    :param str condition: Experimental condition analysed, for details please
    see https://www.ebi.ac.uk/ols
    :param float fpkm_cutoff: threshold for detecting if a gene is expressed or
    not
    :return list: list of filtered genes
    """

    try:
        organism = assembly_to_organism[assembly]
    except KeyError:
        log.error("Unrecognised assembly %s" % assembly)
        raise

    filtered_studies = retrieve.get_run_by_organism(
        organism=organism, condition=condition)
    if condition and not filtered_studies:
        print(read_table({}))
        log.error('No studies found with %s condition' % condition)
        raise ValueError

    filtered_studies = [study for study in filtered_studies if
                        study['ASSEMBLY_USED'] == assembly and
                        study['STATUS'] == 'Complete']
    studies_ids = [study['STUDY_ID'] for study in filtered_studies]

    studies = retrieve.get_studies_by_organism(organism)
    studies = [study for study in studies if
               study['STUDY_ID'] in studies_ids]

    try:
        data = concat([
            read_table(study['GENES_FPKM_COUNTS_FTP_LOCATION'])
            for study in studies])
    except ValueError:
        log.error('No studies found for %s after filtering' % assembly)
        raise
    genes = data.pop('Genes ID')
    return genes[data.max(axis=1) > fpkm_cutoff].unique.tolist()


if __name__ == "__main__":
    pass
