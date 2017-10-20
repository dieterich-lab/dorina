#!/usr/bin/env python
# -*- coding: utf-8
"""
Created on 16:08 11/10/2017 2017 

"""
from subprocess import check_call

import sys

import click
import logging

from dorina import __version__
from dorina.config import config
from dorina.ensembl import EnsemblGFF

log = logging.getLogger('dorina.config')


@click.version_option(version=__version__, prog_name='dorina')
@click.group(context_settings={'help_option_names': ['-h', '--help']})
def cli():
    """
    This is the utility library of the doRiNA database of posttranscriptional regulatory elements.
    """
    pass


@click.command()
@click.option('-r', '--release', type=str, default=config['DEFAULT'].get('version'),
              help='Ensembl release version', show_default=True)
@click.option('-o', '--organism', type=str, default=config['DEFAULT'].get('organism'),
              help='Organism name', show_default=True)
# @click.option('--variants', is_flag=True, help='Retrieves variants')
# @click.option('--dif_expression')
# @click.option('--regulatory')
def create_assembly(release, organism):
    """Retrieves the data files for a given Ensembl release version and organism"""
    log.info('Retrieving gff files for {} {}'.format(release, organism))
    gff = EnsemblGFF(organism=organism, version=release)
    # assembly = gff.assembly
    assembly = 'GRCh38'
    base_path = config['DEFAULT'].get('data_path')
    log.info('Setting up data directory for {} assembly: {}'.format(assembly, base_path))
    # with open("{}/{}.genome".format(base_path, assembly), 'w') as open_f:
    #     check_call(
    #         'mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e'.split() +
    #         ['"select chrom, size from {}.chromInfo"'.format(assembly.replace('GRC', '') )],
    #         stdout=open_f)
    #
    mapping = {"gene": "all",
               "CDS": "cds",
               "exon": "exon",
               "three_prime": "3_utr",
               "five_prime": "5_utr"}

    for k, v in mapping.items():
        with open("{}/{}.bed".format(assembly, k), 'w') as open_f:
            check_call(
                'grep {} {}.gff3'.format(assembly, v).split(),
                stdout=open_f)

    with open("{}/intergenic.bed".format(assembly), 'w') as open_f:
        check_call(
            'bedtools complement -i genes.gff -g {}.genome'.format(assembly).split(),
            stdout=open_f)

    with open("{}/intron.gff3 ".format(assembly), 'w') as open_f:
        check_call(
            'bedtools subtract -s -a genes.gff -b exon.gff |sed -e "s/\tgene\t/\tintron\t/"'.split(),
            stdout=open_f)


cli.add_command(create_assembly)
if __name__ == '__main__':
    cli()
