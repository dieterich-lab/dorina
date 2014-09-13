# vim: set fileencoding=utf-8 :

import os
import logging
from os import path
from cStringIO import StringIO
from pybedtools import BedTool
from dorina import utils


def analyse(genome, set_a, match_a='any', region_a='any',
            set_b=None, match_b='any', region_b='any',
            combine='or', genes=None, slop=0,
            datadir=None):
    """Run doRiNA analysis"""
    logging.debug("analyse(%r, %r(%s) <-'%s'-> %r(%s))" % (genome, set_a, match_a, combine, set_b, match_b))

    return _parse_results(_analyse(genome, set_a, match_a, region_a,
                                   set_b, match_b, region_b, combine,
                                   genes, slop, datadir))


def _analyse(genome, set_a, match_a='any', region_a='any',
             set_b=None, match_b='any', region_b='any', combine='or',
             genes=None, slop=0, datadir=None):
    """Run doRiNA analysis, internal logic"""
    genome_bed_a = _get_genome_bedtool(genome, region_a, datadir, genes)
    regulators_a = map(lambda x: _get_regulator_bedtool(x, datadir), set_a)
    if slop > 0:
        regulators_a = map(lambda x: _add_slop(x, genome, slop, datadir), regulators_a)

    if match_a == 'any':
        regulator = _merge_regulators(regulators_a)
        result_a = _cleanup_intersect_gff(genome_bed_a.intersect(regulator, wa=True, wb=True))
    elif match_a == 'all':
        results_a = map(lambda x: _cleanup_intersect_gff(genome_bed_a.intersect(x, wa=True, wb=True)), regulators_a)
        result_a = None
        for res in results_a:
            if result_a is None:
                result_a = res
                continue
            result_a = _cleanup_intersect_gff_gff(result_a.intersect(res, wa=True, wb=True))

    if set_b is not None:
        genome_bed_b = _get_genome_bedtool(genome, region_b, datadir, genes)
        regulators_b = map(lambda x: _get_regulator_bedtool(x, datadir), set_b)
        if slop > 0:
            regulators_b = map(lambda x: _add_slop(x, genome, slop, datadir), regulators_b)

        if match_b == 'any':
            regulator = _merge_regulators(regulators_b)
            result_b = _cleanup_intersect_gff(genome_bed_b.intersect(regulator, wa=True, wb=True))
        elif match_b == 'all':
            results_b = map(lambda x: _cleanup_intersect_gff(genome_bed_b.intersect(x, wa=True, wb=True)), regulators_b)
            result_b = None
            for res in results_b:
                if result_b is None:
                    result_b = res
                    continue
                result_b = _cleanup_intersect_gff_gff(result_b.intersect(res, wa=True, wb=True))


        if combine == 'or':
            final_results = _merge_regulators([result_a, result_b])
        elif combine == 'and':
            final_results = _cleanup_intersect_gff_gff(result_a.intersect(result_b, wa=True, wb=True))
        elif combine == 'xor':
            not_in_b = result_a.intersect(result_b, v=True, wa=True)
            not_in_a = result_b.intersect(result_a, v=True, wa=True)
            final_results = _merge_regulators([not_in_b, not_in_a])
        elif combine == 'not':
            final_results = result_a.intersect(result_b, v=True, wa=True)

    else:
        final_results = result_a

    return final_results

def _merge_regulators(regulators):
    """Merge a list of regulators using BedTool.cat"""
    regulator = regulators[0]
    for i in range(1, len(regulators)):
        logging.debug('merging regulator %r' % regulators[i])
        regulator = regulator.cat(regulators[i], postmerge=False)

    return regulator


def _intersect_regulators(regulators):
    """Intersect a list of regulators using BedTool.intersect"""
    regulator = regulators[0]
    for i in range(1, len(regulators)):
        logging.debug('intersect regulator %r' % regulators[i])
        dirty_reg = regulator.intersect(regulators[i], wa=True, wb=True)
        regulator = _cleanup_intersect_bed(dirty_reg)

    return regulator


def _cleanup_intersect_bed(dirty):
    clean_string = ''
    for row in dirty:
        try:
            # Bed9 format?
            new_start = "%s" % max(int(row[1]), int(row[9]))
            new_end = "%s" % min(int(row[2]), int(row[10]))
            combined_start = "%s~%s" % (row[6], row[9])
            combined_end = "%s~%s" % (row[7], row[10])
            new_name = "~".join([row[3], row[11]])
            new_score = "%s~%s" % (row[4], row[12])
            new_strand = row[5] if row[5] == row[13] else '.'
        except ValueError:
            # Try bed6 instead
            new_start = "%s" % max(int(row[1]), int(row[7]))
            new_end = "%s" % min(int(row[2]), int(row[8]))
            combined_start = "%s~%s" % (row[1], row[7])
            combined_end = "%s~%s" % (row[2], row[8])
            new_name = "~".join([row[3], row[9]])
            new_score = "%s~%s" % (row[4], row[10])
            new_strand = row[5] if row[5] == row[11] else '.'

        new_row = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(
            row[0], new_start, new_end, new_name, new_score, new_strand,
            combined_start, combined_end)
        clean_string += new_row

    return BedTool(clean_string, from_string=True)


def _cleanup_intersect_gff(dirty):
    clean_string = ''
    for row in dirty:
        new_row = "\t".join(row[:8])
        annotations = _parse_annotations(row[8])

        if len(row.fields) >= 17:
            if 'start' in annotations:
                annotations['start'] = "{0}~{1}".format(annotations['start'], row[15])
                annotations['end'] = "{0}~{1}".format(annotations['end'], row[16])
            else:
                annotations['start'] = row[15]
                annotations['end'] = row[16]
        else:
            if 'start' in annotations:
                annotations['start'] = "{0}~{1}".format(annotations['start'], row[10])
                annotations['end'] = "{0}~{1}".format(annotations['end'], row[11])
            else:
                annotations['start'] = row[10]
                annotations['end'] = row[11]

        if 'regulator' in annotations:
            annotations['regulator'] = "{0}~{1}".format(annotations['regulator'], row[12])
        else:
            annotations['regulator'] = row[12]

        if 'score' in annotations:
            annotations['score'] = "{0}~{1}".format(annotations['score'], row[13])
        else:
            annotations['score'] = row[13]

        annotations['ID'] += "~{}".format(annotations['ID'])

        new_row += "\tID={ID};regulator={regulator};score={score};start={start};end={end}\n".format(**annotations)
        clean_string += new_row

    return BedTool(clean_string, from_string=True)

def _cleanup_intersect_gff_gff(dirty):
    clean_string = ''
    for row in dirty:
        new_row = "\t".join(row[:8])
        ann_a = _parse_annotations(row[8])
        ann_b = _parse_annotations(row[17])
        new_annotations = {}
        new_annotations['regulator'] = "{0}~{1}".format(ann_a['regulator'], ann_b['regulator'])
        new_annotations['score'] = "%s~%s" % (ann_a['score'], ann_b['score'])
        new_annotations['start'] = "%s~%s" % (ann_a['start'], ann_b['start'])
        new_annotations['end'] = "%s~%s" % (ann_a['end'], ann_b['end'])
        new_annotations['ID'] = "{0}~{1}".format(ann_a['ID'], ann_b['ID'])
        new_row += "\tID={0};regulator={1};score={2};start={3};end={4}\n".format(new_annotations['ID'],
            new_annotations['regulator'], new_annotations['score'],
            new_annotations['start'], new_annotations['end'])
        clean_string += new_row

    return BedTool(clean_string, from_string=True)


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
        res = filter_name in rec.name
        return res

    if os.sep in regulator_name:
        return BedTool('%s.bed' % utils.get_regulator_by_name(regulator_name, datadir))

    return BedTool('%s.bed' % utils.get_regulator_by_name(regulator_name, datadir)).filter(filter_func, regulator_name).saveas()


def _parse_results(bedtool_results):
    """parse a bedtool result data structure"""
    results = []

    for res in bedtool_results:
        annotations = _parse_annotations(res[8])
        genes = annotations['ID'].split('~')
        tracks, data_sources, sites = _parse_tracks_sources_regulators(annotations['regulator'])
        scores = annotations['score'].split('~')
        starts = annotations['start'].split('~')
        ends = annotations['end'].split('~')
        for i in range(len(tracks)):
            gene = genes[i]
            track = tracks[i]
            data_source = data_sources[i]
            site = sites[i]
            start = starts[i]
            end = ends[i]
            strand = res.strand
            try:
                score = float(scores[i])
            except ValueError:
                logging.debug("can't parse score '{}'".format(scores[i]))
                score = 0

            location = "%s:%s-%s" % (res.chrom, start, end)
            new_res = dict(track=track, gene=gene, data_source=data_source,
                           score=score, site=site, location=location, strand=strand)
            if new_res not in results:
                results.append(new_res)

    return results


def _parse_annotations(string):
    annotation_dict = {}
    annotation_list = string.split(';')
    for ann in annotation_list:
        key, val = ann.split('=')
        annotation_dict[key] = val

    return annotation_dict


def _parse_tracks_sources_regulators(string):
    raw = string.split('~')
    tracks = []
    sources = []
    regulators = []
    for r in raw:
        try:
            source, rest = r.split('#')
        except ValueError:
            source = "Unknown"
            rest = r

        try:
            track, regulator = rest.split('*')
        except ValueError:
            track = "Unknown"
            regulator = rest

        tracks.append(track)
        sources.append(source)
        regulators.append(regulator)

    return tracks, sources, regulators
