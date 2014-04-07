# vim: set fileencoding=utf-8 :

import logging
from os import path
from cStringIO import StringIO
from pybedtools import BedTool
from dorina import utils


def analyse(genome, set_a, match_a='any', region_a='any', genes=None, datadir=None):
    """Run doRiNA analysis"""
    logging.debug("analyse(%r, %r(%s))" % (genome, set_a, match_a))

    return _parse_results(_analyse(genome, set_a, match_a, region_a, datadir))


def _analyse(genome, set_a, match_a='any', region_a='any', datadir=None):
    """Run doRiNA analysis, internal logic"""
    logging.debug("analyse(%r, %r(%s))" % (genome, set_a, match_a))

    genome_bed = _get_genome_bedtool(genome, region_a, datadir)
    regulators = map(lambda x: _get_regulator_bedtool(x, datadir), set_a)

    if match_a == 'any':
        regulator = _merge_regulators(regulators)
    elif match_a == 'all':
        regulator = _intersect_regulatrs(regulators)

    return genome_bed.intersect(regulator, wa=True, wb=True)


def _merge_regulators(regulators):
    """Merge a list of regulators using BedTool.cat"""
    regulator = regulators[0]
    for i in range(1, len(regulators)):
        logging.debug('merging regulator %r' % regulators[i])
        regulator = regulator.cat(regulators[i], postmerge=False)

    return regulator


def _intersect_regulatrs(regulators):
    """Intersect a list of regulators using BedTool.intersect"""
    regulator = regulators[0]
    for i in range(1, len(regulators)):
        logging.debug('intersect regulator %r' % regulators[i])
        regulator = regulator.intersect(regulators[i], wa=True, wb=True)

    return regulator


def _get_genome_bedtool(genome_name, region, datadir=None):
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

    return BedTool(filename)


def _get_regulator_bedtool(regulator_name, datadir=None):
    """get the bedtool object for a regulator"""
    return BedTool('%s.bed' % utils.get_regulator_by_name(regulator_name, datadir))


def _parse_results(bedtool_results):
    """parse a bedtool result data structure"""
    results = []

    for res in bedtool_results:
        line = str(res).strip().split('\t')

        track = line[0]
        gene = line[8]
        data_source = line[12].split('_')[0]
        score = int(line[13])
        site = line[12]
        strand = line[6]
        start = line[10]
        end = line[11]
        if len(line) > 19:
            start = min(start, line[18])
            end = max(end, line[19])
        location = "%s:%s-%s" % (track, start, end)

        results.append(dict(track=track, gene=gene, data_source=data_source,
                            score=score, site=site, location=location, strand=strand))

    return results
