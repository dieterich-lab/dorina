# vim: set fileencoding=utf-8 :

import os
import logging
from os import path
from pybedtools import BedTool
from dorina import utils


def analyse(genome, set_a, match_a='any', region_a='any',
            set_b=None, match_b='any', region_b='any',
            combine='or', genes=None, slop=0,
            datadir=None):
    """Run doRiNA analysis"""
    logging.debug("analyse(%r, %r(%s) <-'%s'-> %r(%s))" % (genome, set_a, match_a, combine, set_b, match_b))

    def compute_result(region, regulators, match):
        genome_bed = _get_genome_bedtool(genome, region, datadir, genes)
        if slop > 0:
            regulators = map(lambda x: _add_slop(x, genome, slop, datadir), regulators)

        result = None
        if match == 'any':
            result = genome_bed.intersect(_merge_regulators(regulators), wa=True, u=True)
        elif match == 'all':
            results = [genome_bed]
            results.extend(regulators)
            result = reduce(lambda acc, x: acc.intersect(x, wa=True, u=True), results)

        return result

    regulators_a = map(lambda x: _get_regulator_bedtool(x, datadir), set_a)
    regulators_b = []

    result_a = compute_result(region_a, regulators_a, match_a)
    if set_b is not None:
        regulators_b = map(lambda x: _get_regulator_bedtool(x, datadir), set_b)
        result_b = compute_result(region_b, regulators_b, match_b)
        if combine == 'or':
            final_results = _merge_regulators([result_a, result_b])
        elif combine == 'and':
            final_results = result_a.intersect(result_b, wa=True, u=True)
        elif combine == 'xor':
            not_in_b = result_a.intersect(result_b, v=True, wa=True)
            not_in_a = result_b.intersect(result_a, v=True, wa=True)
            final_results = _merge_regulators([not_in_b, not_in_a])
        elif combine == 'not':
            final_results = result_a.intersect(result_b, v=True, wa=True)

    else:
        final_results = result_a

    regulators_a.extend(regulators_b)
    regulator = _merge_regulators(regulators_a)
    final_results = final_results.intersect(regulator, wa=True, wb=True)
    return final_results

def _merge_regulators(regulators):
    """Merge a list of regulators using BedTool.cat"""
    regulator = regulators[0]
    for i in range(1, len(regulators)):
        logging.debug('merging regulator %r' % regulators[i])
        regulator = regulator.cat(regulators[i], postmerge=False)

    return regulator


def _add_slop(feature, genome_name, slop, datadir=None):
    """Add specified slop before and after a regulator"""
    return feature.slop(g=_get_genome_chromfile(genome_name, datadir), b=slop)


def _get_genome_chromfile(genome_name, datadir=None):
    """Get the path to the .genome file listing chromosome sizes"""
    return path.join(utils.get_genome_by_name(genome_name, datadir),
                     "{}.genome".format(genome_name))


def _get_genome_bedtool(genome_name, region, datadir=None, genes=None):
    """get the bedtool object for a genome depending on the name and the region"""
    genome = utils.get_genome_by_name(genome_name, datadir)
    if region == "any":
        filename = path.join(genome, 'all.gff')
    elif region == "CDS":
        filename = path.join(genome, 'cds.gff')
    elif region == "3prime":
        filename = path.join(genome, '3_utr.gff')
    elif region == "5prime":
        filename = path.join(genome, '5_utr.gff')
    elif region == "intron":
        filename = path.join(genome, 'intron.gff')
    elif region == "intergenic":
        filename = path.join(genome, 'intergenic.gff')
    else:
        raise ValueError("Invalid region: %r" % region)

    if genes is None or 'all' in genes:
        return BedTool(filename)
    else:
        return BedTool(filename).filter(lambda x: x.name in genes).saveas()

def _get_regulator_bedtool(regulator_name, datadir=None):
    """get the bedtool object for a regulator"""
    def filter_func(rec, name):
        if '_' in name:
            filter_name = '_'.join(name.split('_')[1:])
        else:
            filter_name = name
        res = (filter_name + '*' in rec.name) or (filter_name == rec.name)
        return res

    if os.sep in regulator_name:
        return BedTool('%s.bed' % utils.get_regulator_by_name(regulator_name, datadir))

    return BedTool('%s.bed' % utils.get_regulator_by_name(regulator_name, datadir)).filter(filter_func, regulator_name).saveas()
