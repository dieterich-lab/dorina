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


            result = None
            if match == 'any':
                result = genome_bed.intersect(self._merge_regulators(_regulators), wa=True, u=True)
            elif match == 'all':
                results = [genome_bed]
                results.extend(_regulators)
                result = reduce(lambda acc, x: acc.intersect(x, wa=True, u=True), results)

            return result

        regulators_a = map(lambda x: self._get_regulator_bedtool(x), set_a)
        regulators_b = []

        result_a = compute_result(region_a, regulators_a, match_a, window_a)
        if set_b is not None:
            regulators_b = map(lambda x: self._get_regulator_bedtool(x), set_b)
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
        return feature.slop(g=self._get_genome_chromfile(genome_name), b=slop)

    def _get_genome_chromfile(self, genome_name):
        """Get the path to the .genome file listing chromosome sizes"""
        return path.join(self.utils.get_genome_by_name(genome_name),
                         "{}.genome".format(genome_name))


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
            filename = path.join(genome, "%s.gff" % mapping[region])

        if genes is None or 'all' in genes:
            return BedTool(filename)
        else:
            return BedTool(filename).filter(lambda x: x.name in genes).saveas()

    def _get_regulator_bedtool(self, regulator_name):
        """get the bedtool object for a regulator"""
        def filter_func(rec, name):
            if '_' in name:
                filter_name = '_'.join(name.split('_')[1:])
            else:
                filter_name = name
            res = (filter_name + '*' in rec.name) or (filter_name == rec.name)
            return res

        if os.sep in regulator_name or '_all' in regulator_name:
            bt = BedTool('%s.bed' % self.utils.get_regulator_by_name(regulator_name))
        else:
            bt = BedTool('%s.bed' % self.utils.get_regulator_by_name(regulator_name)).filter(filter_func, regulator_name).saveas()

        if len(bt) > 0 and len(bt[0].fields) > 6:
            bt = bt.bed6().saveas()

        return bt
