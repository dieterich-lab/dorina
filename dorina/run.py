# vim: set fileencoding=utf-8 :

import os
import logging
from os import path
from pybedtools import BedTool

from dorina import utils

class Dorina:
    def __init__(self, datadir):
        self.utils = utils.DorinaUtils(datadir)

    def analyse(self, genome,
                set_a,      match_a='any', region_a='any',
                set_b=None, match_b='any', region_b='any',
                combine='or', genes=None,
                window_a=-1,
                window_b=-1):
        """Run doRiNA analysis"""
        logging.debug("analyse(%r, %r(%s) <-'%s'-> %r(%s))" % (genome, set_a, match_a, combine, set_b, match_b))

        def compute_result(region, regulators, match, window):
            genome_bed = self._get_genome_bedtool(genome, region, genes)

            # create local copy so we can mangle it
            _regulators = regulators[:]
            if window > -1:
                initial = _regulators.pop(0)
                genome_bed = genome_bed.intersect(initial)
                if window > 0:
                    genome_bed = self._add_slop(genome_bed, genome, window)

            if match == 'any':
                result = genome_bed.intersect(self._merge_regulators(_regulators), wa=True, u=True)
            elif match == 'all':
                result = reduce(lambda acc, x: acc.intersect(x, wa=True, u=True),
                                [genome_bed] + _regulators)
            else:
                result = None
            return result

        regulators_a = self._regulators_from_names(set_a)
        regulators_b = self._regulators_from_names(set_b)

        result_a = compute_result(region_a, regulators_a, match_a, window_a)

        # Combine with set B, if exists
        if set_b:
            result_b = compute_result(region_b, regulators_b, match_b, window_b)
            if combine == 'or':
                final_results = self._merge_regulators([result_a, result_b])
            elif combine == 'and':
                final_results = result_a.intersect(result_b, wa=True, u=True)
            elif combine == 'xor':
                not_in_b = result_a.intersect(result_b, v=True, wa=True)
                not_in_a = result_b.intersect(result_a, v=True, wa=True)
                final_results = self._merge_regulators([not_in_b, not_in_a])
            elif combine == 'not':
                final_results = result_a.intersect(result_b, v=True, wa=True)

        else:
            final_results = result_a

        regulators_a.extend(regulators_b)
        regulator = self._merge_regulators(regulators_a)
        final_results = final_results.intersect(regulator, wa=True, wb=True)
        return final_results

    def _merge_regulators(self, regulators):
        """Merge a list of regulators using BedTool.cat"""
        if len(regulators) > 1:
            return BedTool.cat(*regulators, postmerge=False)
        else:
            return regulators[0]

    def _add_slop(self, feature, genome_name, slop):
        """Add specified slop before and after a regulator"""
        genome = self.utils.get_genome_by_name(genome_name)
        chromfile = path.join(genome, "{}.genome".format(genome_name))
        return feature.slop(g=chromfile, b=slop)

    def _get_genome_bedtool(self, genome_name, region, genes=None):
        """get the bedtool object for a genome depending on the name and the region"""
        genome = self.utils.get_genome_by_name(genome_name)
        mapping = { "any":        "all",
                    "CDS":        "cds",
                    "3prime":     "3_utr",
                    "5prime":     "5_utr",
                    "intron":     "intron",
                    "intergenic": "intergenic" }

        if region not in mapping:
            raise ValueError("Invalid region: %r" % region)
        else:
            bed = BedTool(path.join(genome, "%s.gff" % mapping[region]))

        # Optionally, filter by gene.
        if genes is None or 'all' in genes:
            return bed
        else:
            return bed.filter(lambda x: x.name in genes).saveas()

    # TODO: write regulator class
    def _regulators_from_names(self, names):
        if names:
            return map(lambda x: self._get_regulator_bedtool(x), names)
        else:
            return []

    def _get_regulator_bedtool(self, regulator_name):
        """get the bedtool object for a regulator"""

        def by_name(rec):
            name = regulator_name
            # Drop first part before underscore.
            if "_" in name:
                name = "_".join(name.split("_")[1:])
            return (name + "*" in rec.name) or (name == rec.name)

        bt = BedTool('%s.bed' % self.utils.get_regulator_by_name(regulator_name))
        if os.sep not in regulator_name and '_all' not in regulator_name:
            bt = bt.filter(by_name).saveas()

        if len(bt) > 0 and len(bt[0].fields) > 6:
            bt = bt.bed6().saveas()

        return bt
