#!/usr/bin/env python
# -*- coding: utf-8
"""
Created on 16:08 11/10/2017 2017
"""
# pragma: no cover
from __future__ import unicode_literals

import functools
import logging
import os
import shutil
import sys
from io import open
from subprocess import check_call, Popen, PIPE

import click

from dorina import __version__, run as run_dorina
from dorina.config import config
from dorina.ensembl import EnsemblFTP
from dorina.genome import Genome
from dorina.regulator import Regulator

# https://stackoverflow.com/a/15729700/1694714
log = logging.getLogger(__name__)
logging.captureWarnings(True)
logging.basicConfig(stream=sys.stderr, level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s ')


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
    database of posttranscriptional regulatory elements
    """
    pass


def common_params(func):
    @click.option('-r', '--release', type=str,
                  default=config.get('DEFAULT', 'version'),
                  help='Ensembl release version', show_default=True)
    @click.option('-o', '--organism',
                  type=str, default=config.get('DEFAULT', 'organism'),
                  help='Organism formal name', show_default=True)
    @click.option('--variation', is_flag=True,
                  help='Retrieves data set from Ensembl variation')
    @click.option('--regulation', is_flag=True,
                  help='Retrieves data set from Ensembl regulation')
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


@click.command()
@common_params
def create_assembly(release, organism, variation, regulation):
    """Retrieves the data files for a given Ensembl release version and
    organism """
    log.info('Retrieving gff files for {} {}'.format(release, organism))
    ftp = EnsemblFTP(organism=organism, release=release)
    try:
        assembly = ftp.assembly
    except (IndexError,):
        raise ValueError('No assembly found for {}.'.format(
            organism))
    base_path = config.get('DEFAULT', 'data_path')

    gff = ftp.retrieve_full_gff()[0]

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

    for k, v in list(mapping.items()):
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
        "bedtools complement -i all_sorted.gff3 -g {}.genome".format(
            assembly).split(),
        stdout="intergenic.bed",
        cwd=os.path.join(base_path, assembly))

    os.remove(os.path.join(base_path, assembly, 'exon.gff3'))

    if variation:
        ftp.retrieve_full_vcf()

    if regulation:
        ftp.retrieve_from_regulation_by_experiment()

    sys.exit(0)


@click.command()
@click.argument('assembly', type=str)
def clear_assembly(assembly):
    """Clear data files for a given Ensembl release version and organism """
    click.confirm(
        'Do you want to remove data file for the {} assembly?'.format(assembly),
        abort=True)

    base_path = config.get('DEFAULT', 'data_path')
    shutil.rmtree(os.path.join(base_path, assembly))


@click.command()
@click.argument('genome')
@click.option('-d', '--debug', is_flag=True,
              help="Set logging level to debug (more verbose)")
@click.option('-q', '--quiet', is_flag=True,
              help="Set logging level to error (quieter)")
@click.option('-a', '--seta', required=True, multiple=True,
              help="First set of regulators to analyse")
@click.option('-b', '--setb',
              help="Second set of regulators to analyse")
@click.option('--genes', multiple=True, default=['all'])
@click.option('--matcha', required=True, type=click.Choice(['any', 'all']),
              default='any', help="All or any regulators in set A must match")
@click.option('--regiona', default='any', type=click.Choice(
    ['any', 'CDS', '3prime', '5prime', 'intron', 'intergenic']),
              help="Region to match set A in")
@click.option('--matchb', type=click.Choice(['any', 'all']), default='any',
              help="All or any regulators in set B must match")
@click.option('--regionb', default='any', type=click.Choice(
    ['any', 'CDS', '3prime', '5prime', 'intron', 'intergenic']),
              help="Region to match set B in")
@click.option('-C', '--combine', default='or',
              type=click.Choice(['and', 'or', 'not', 'xor']),
              help="Set operation to combine set A and set B hits")
@click.option('--windowa', type=int, default=-1,
              help="Use windowed search for set A")
@click.option('--windowb', type=int, default=-1,
              help="Use windowed search for set B")
@click.option('--path', '-p', default=config.get('DEFAULT', 'data_path'),
              type=click.Path(exists=True, dir_okay=True, readable=True),
              help="Path to genomes and regulators")
def run(genome, debug, quiet, seta, setb, genes, matcha, regiona,
        matchb, regionb, combine, windowa, windowb, path):
    """"Run doRiNA from the command line"""
    if debug:
        log.setLevel(logging.DEBUG)
    elif quiet:
        log.setLevel(logging.ERROR)

    dorina = run_dorina.Dorina(path)
    Genome.init(path)
    mapping = {}
    for x in Genome.all().values():
        for y in x['assemblies']:
            mapping[y] = x['id']

    # if 'all' in seta:  # seta is a tuple
    #     seta = Regulator.all()[mapping[genome]][genome].keys()

    result = dorina.analyse(genome, seta, matcha, regiona, setb, matchb,
                            regionb, combine, genes, windowa, windowb)
    click.echo(result)
    sys.exit(0)


@click.command()
@click.option('--path', '-p', default=config.get('DEFAULT', 'data_path'),
              type=click.Path(exists=True, dir_okay=True, readable=True),
              help="Path to genomes and regulators")
def genomes(path):
    """List available genomes in given directory"""
    Genome.init(path)
    _genomes = Genome

    if _genomes is None:
        click.echo('No genomes available.')
        sys.exit()

    click.echo("Available genomes:")
    click.echo("------------------")
    for species, species_dict in _genomes.all().items():
        click.echo("\t%s" % species)
        for assembly, assembly_dict in species_dict['assemblies'].items():
            click.echo("\t\t%s" % assembly)
            gffs = list(assembly_dict.items())
            gffs.sort(key=lambda x: x[0])
            for gff in gffs:
                click.echo("\t\t\t%s: %s" % gff)
    sys.exit(0)


@click.command()
@click.option('--path', '-p', default=config.get('DEFAULT', 'data_path'),
              type=click.Path(exists=True, dir_okay=True, readable=True),
              help="Path to genomes and regulators")
def regulators(path):
    """List available regulators in a given directory"""
    Regulator.init(path)
    _regulator = Regulator.all()
    if _regulator is None:
        click.echo('No regulators available.')
        sys.exit()

    click.echo("Available regulators:")
    click.echo("---------------------")
    for species, species_dict in _regulator.items():
        click.echo("\t%s" % species)
        for assembly, assembly_dict in species_dict.items():
            click.echo("\t\t%s" % assembly)
            for _regulator, regulator_dict in list(assembly_dict.items()):
                click.echo("\t\t\t%s" % _regulator)
    sys.exit(0)


cli.add_command(create_assembly)
cli.add_command(clear_assembly)
cli.add_command(regulators)
cli.add_command(genomes)
cli.add_command(run)
if __name__ == '__main__':
    cli("run hg19 -a all -p /Volumes/prj/dorina2/".split())
