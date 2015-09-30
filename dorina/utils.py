# vim: set fileencoding=utf-8 :

import logging
import os
import json
from os import path
import re
from pybedtools import BedTool

gene_name = re.compile(r'.*ID=(.*?)($|;\w+)')

class DorinaUtils:
    def __init__(self, datadir):
        self.datadir = datadir
        self.genomes = self._genomes()
        self.regulators = self._regulators()

    def _genomes(self):
        """Get all available genomes"""
        def parse_func(assembly_path, assembly_dict):
            for gff_file in os.listdir(assembly_path):
                gff_path = path.join(assembly_path, gff_file)
                if not path.isfile(gff_path):
                    #logging.debug("skipping non-file %r" % gff_path)
                    continue

                root, ext = path.splitext(gff_file)
                if ext in ('.gff', '.bed'):
                    assembly_dict[root] = True

        return self.walk_assembly_tree('genomes', parse_func)

    def _regulators(self):
        """Get all available regulators.  A valid regulator must have a JSON metadata
file as well as a BED file containing the data."""
        def parse_func(root, regulators):
            for experiment in os.listdir(root):
                experiment_path = path.join(root, experiment)
                if not path.isfile(experiment_path):
                    #logging.debug("skipping non-file %r" % experiment_path)
                    continue

                experiment_root, experiment_ext = path.splitext(experiment)
                if not experiment_ext.lower() == '.json':
                    #logging.debug("skipping non-JSON file %r" % experiment_path)
                    continue

                bedfile = path.join(root, '%s.%s' % (experiment_root, 'bed'))
                #logging.debug("looking for %r" % bedfile)
                if not path.isfile(bedfile):
                    #logging.debug("No bedfile for experiment %r" % experiment_path)
                    continue

                experiments = self.parse_experiment(experiment_path)
                for experiment_dict in experiments:
                    experiment_dict['file'] = experiment_path
                    regulators[experiment_dict['id']] = experiment_dict

        return self.walk_assembly_tree('regulators', parse_func)

    def parse_experiment(self, filename):
        """Parse experimental description from a file name"""
        #logging.debug("Parsing experimental description from %r" % filename)
        experiment = {}
        with open(filename, 'r') as fh:
            experiment = json.load(fh)

        return experiment

    def walk_assembly_tree(self, root_dir, parse_func):
        """Walk a directory structure containg clade, species, assembly

        Call parse_func() for every assembly directory"""
        genomes = {}

        root = path.join(self.datadir, root_dir)

        for species in os.listdir(root):
            species_path = path.join(root, species)
            if not path.isdir(species_path):
                #logging.debug("skipping non-directory %r" % species_path)
                continue

            description_file = path.join(species_path, 'description.json')
            species_dict = {}
            try:
                with open(description_file, 'r') as fh:
                    genomes[species] = json.load(fh)
                    genomes[species]['assemblies'] = species_dict
            except IOError:
                genomes[species] = species_dict

            for assembly in os.listdir(species_path):
                assembly_path = path.join(species_path, assembly)
                if not path.isdir(assembly_path):
                    #logging.debug("skipping non-directory %r" % assembly_path)
                    continue

                assembly_dict = {}
                species_dict[assembly] = assembly_dict

                parse_func(assembly_path, assembly_dict)

        return genomes

    def get_genome_by_name(self, name):
        """Take a genome name and return the path to the genome directory"""
        for species, species_dir in self.genomes.items():
            if name in species_dir['assemblies']:
                return path.join(self.datadir, 'genomes', species, name)

        # TODO: never return None!
        return None

    def get_genes(self, name):
        """Get a list of genes from genome <name>"""
        genes = []

        genome_dir = self.get_genome_by_name(name)
        if genome_dir is None:
            return genes

        genome = path.join(genome_dir, 'all.gff')
        if not path.exists(genome):
            return genes

        for line in open(genome, 'r'):
            match = gene_name.match(line)
            if match is not None:
                genes.append(match.group(1))

        return genes

    def make_regulator(self, name_or_path, assembly=None):
        if os.sep in name_or_path:
            return Regulator("custom", name_or_path, True)

        if not assembly:
            raise ValueError("Must provide assembly")

        filename = None
        for species, species_dir in self.regulators.items():
            for _assembly, assembly_dir in species_dir.items():
                if assembly == _assembly and name_or_path in assembly_dir:
                    basename = path.splitext(assembly_dir[name_or_path]['file'])[0]
                    filename = basename + ".bed"

        if not filename or not path.isfile(filename):
            raise ValueError("Could not find regulator: %s" % name_or_path)

        return Regulator(name_or_path, filename, False)

    def regulators_from_names(self, names, assembly):
        if names:
            return map(lambda x: self.make_regulator(x, assembly).bed, names)
        else:
            return []


class Regulator:
    def __init__(self, name, path, custom):
        self.name = name
        self.path = path
        self.basename = os.path.splitext(path)[0]
        self.custom = custom
        self.bed = self._bed()

    def _bed(self):
        def by_name(rec):
            # Drop first part before underscore.
            if "_" in self.name:
                name = "_".join(self.name.split("_")[1:])
            else:
                name = self.name
            return (name + "*" in rec.name) or (name == rec.name)

        bt = BedTool(self.path)
        if not self.custom and '_all' not in self.name:
            bt = bt.filter(by_name).saveas()

        if len(bt) > 0 and len(bt[0].fields) > 6:
            bt = bt.bed6().saveas()

        return bt
