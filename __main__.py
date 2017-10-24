#!/usr/bin/env python
# -*- coding: utf-8
"""
Created on 16:08 11/10/2017 2017 

"""
import functools
import logging
import os
from subprocess import check_call, Popen, PIPE

import click
import sys

import shutil

from dorina import __version__
from dorina.config import config
from dorina.ensembl import EnsemblFTP

log = logging.getLogger('dorina.config')


def call_command(command, stdout=None, cwd=None, mode='w', stdin=None):
    """
    Pipe the output of command to basename in the cwd.

    :param stdin:
    :param str stdout: basename for the output file in cwd.
    :param list command: list of executable plus arguments
    :param cwd: current working directory
    :param mode: mode for file output
    """
    if not isinstance(stdout, str):
        return Popen(command, cwd=cwd, stdout=stdout, stdin=stdin)

    target = os.path.join(cwd, stdout)
    try:
        with open(target, mode=mode) as open_f:
            check_call(command, cwd=cwd, stdout=open_f, stdin=stdin)
    except Exception as e:
        os.remove(target)
        raise e


@click.version_option(version=__version__, prog_name='dorina')
@click.group(context_settings={'help_option_names': ['-h', '--help']})
def cli():
    """
    This is the command line interface for the doRiNA, a
    database of posttranscriptional regulatory elements.
    """
    pass


def common_params(func):
    @click.option('-r', '--release', type=str,
                  default=config['DEFAULT'].get('version'),
                  help='Ensembl release version', show_default=True)
    @click.option('-o', '--organism',
                  type=str, default=config['DEFAULT'].get('organism'),
                  help='Organism name', show_default=True)
    # @click.option('--variants', is_flag=True, help='Retrieves variants')
    # @click.option('--dif_expression')
    # @click.option('--regulatory')
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


@click.command()
@common_params
def create_assembly(release, organism):
    """Retrieves the data files for a given Ensembl release version and
    organism """
    #  TODO better error if called with mouse / rat or support
    log.info('Retrieving gff files for {} {}'.format(release, organism))
    ftp = EnsemblFTP(organism=organism, version=release)
    assembly = ftp.assembly
    base_path = config['DEFAULT'].get('data_path')

    gff = ftp.retrieve_from_gff_by_index()[0]

    cl = r'$1 == "##sequence-region"{print $2 "\t" $4 - $3 + 1} !/^ *#/{exit;}'
    log.debug('Creating {}.genome file'.format(assembly))
    call_command(
        ['awk', cl, gff.replace('.gz', '')],
        "{}.genome".format(assembly),
        cwd=os.path.join(base_path, assembly))

    log.info('Data on {}.'.format(os.path.join(base_path, assembly)))

    mapping = {"gene": "all",
               "CDS": "cds",
               "exon": "exon",
               "three_prime": "3_utr",
               "five_prime": "5_utr"}

    for k, v in mapping.items():
        log.debug('Creating {}.gff3 file'.format(v))
        call_command(
            "grep {} {}".format(k, gff.replace('.gz', '')).split(),
            stdout=v + ".gff3",
            cwd=os.path.join(base_path, assembly))

    log.debug('Creating intron.gff3 file')
    bedtools_sub = call_command(
        'bedtools subtract -s -a all.gff3 -b exon.gff3'.split(),
        stdout=PIPE,
        cwd=os.path.join(base_path, assembly))
    call_command(
        r"sed -e s/\tgene\t/\tintron\t/".split(' '),  # don't split on \t
        stdout="intron.gff3",
        stdin=bedtools_sub.stdout,
        cwd=os.path.join(base_path, assembly))


    call_command(
        "bedtools sort -i all.gff3".split(),
        stdout="all_sorted.gff3",
        cwd=os.path.join(base_path, assembly))
    os.remove(os.path.join(base_path, assembly, 'all.gff3'))

    log.debug('Creating intergenic.bed file')
    call_command(
        "bedtools complement -i all_sorted.gff3 -g {}.genome".format(assembly).split(),
        stdout="intergenic.bed",
        cwd=os.path.join(base_path, assembly))

    os.remove(os.path.join(base_path, assembly, 'exon.gff3'))
    sys.exit(0)


@click.command()
@click.argument('assembly', type=str)
def clear_assembly(assembly):
    """Clear data files for a given Ensembl release version and organism """
    click.confirm(
        'Do you want to remove data file for the {} assembly?'.format(assembly),
                  abort=True)

    base_path = config['DEFAULT'].get('data_path')
    shutil.rmtree(os.path.join(base_path, assembly))


cli.add_command(create_assembly)
cli.add_command(clear_assembly)
if __name__ == '__main__':
    cli()
