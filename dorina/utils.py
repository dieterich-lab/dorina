# vim: set fileencoding=utf-8 :

import logging
import os
import json
from os import path

def get_genomes(datadir):
    """Get all available genomes"""
    def parse_func(assembly_path, assembly_dict):
        for gff_file in os.listdir(assembly_path):
            gff_path = path.join(assembly_path, gff_file)
            if not path.isfile(gff_path):
                logging.debug("skipping non-file %r" % gff_path)
                continue

            root, ext = path.splitext(gff_file)
            if ext in ('.gff', '.bed'):
                assembly_dict[root] = True

    return walk_assembly_tree(datadir, 'genomes', parse_func)


def get_regulators(datadir):
    """Get all available regulators"""
    def parse_func(root, regulators):
        for regulator in os.listdir(root):
            regulator_path = path.join(root, regulator)
            if not path.isdir(regulator_path):
                logging.debug("skipping non-dir %r" % regulator_path)
                continue

            regulator_dict = {}
            regulators[regulator] = regulator_dict

            for experiment in os.listdir(regulator_path):
                experiment_path = path.join(regulator_path, experiment)
                if not path.isfile(experiment_path):
                    logging.debug("skipping non-file %r" % experiment_path)
                    continue

                experiment_root, experiment_ext = path.splitext(experiment)
                if not experiment_ext.lower() == '.json':
                    logging.debug("skipping non-JSON file %r" % experiment_path)
                    continue

                bedfile = path.join(regulator_path, '%s.%s' % (experiment_root, 'bed'))
                logging.debug("looking for %r" % bedfile)
                if not path.isfile(bedfile):
                    logging.debug("No bedfile for experiment %r" % experiment_path)
                    continue

                experiments = parse_experiment(experiment_path)
                for experiment_dict in experiments:
                    experiment_dict['file'] = experiment_path
                    regulator_dict[experiment_dict['id']] = experiment_dict

    return walk_assembly_tree(datadir, 'regulators',parse_func)


def parse_experiment(filename):
    """Parse experimental description from a file name"""
    logging.debug("Parsing experimental description from %r" % filename)
    experiment = {}
    with open(filename, 'r') as fh:
        experiment = json.load(fh)

    return experiment

def walk_assembly_tree(datadir, root_dir, parse_func):
    """Wa;k a directory structure containg clade, species, assembly

    Call parse_func() for every assembly directory"""
    genomes = {}

    root = path.join(datadir, root_dir)

    for clade in os.listdir(root):
        clade_path = path.join(root, clade)
        if not path.isdir(clade_path):
            logging.debug("skipping non-directory %r" % clade_path)
            continue
        clade_dict = {}
        genomes[clade] = clade_dict

        for species in os.listdir(clade_path):
            species_path = path.join(clade_path, species)
            if not path.isdir(species_path):
                logging.debug("skipping non-directory %r" % species_path)
                continue

            species_dict = {}
            clade_dict[species] = species_dict

            for assembly in os.listdir(species_path):
                assembly_path = path.join(species_path, assembly)
                if not path.isdir(assembly_path):
                    logging.debug("skipping non-directory %r" % assembly_path)
                    continue

                assembly_dict = {}
                species_dict[assembly] = assembly_dict

                parse_func(assembly_path, assembly_dict)

    return genomes


def get_genome_by_name(name, datadir):
    """Take a genome name and return the path to the genome directory"""
    genomes = get_genomes(datadir=datadir)

    for clade, clade_dir in genomes.items():
        for species, species_dir in clade_dir.items():
            if name in species_dir:
                return path.join(datadir, 'genomes', clade, species, name)

    return None


def get_regulator_by_name(name, datadir):
    """Take a regulator name and return the path to the regulator basename without file extension"""
    regulators = get_regulators(datadir=datadir)

    for clade, clade_dir in regulators.items():
        for species, species_dir in clade_dir.items():
            for assembly, assembly_dir in species_dir.items():
                for regulator, regulator_dir in assembly_dir.items():
                    if name in regulator_dir:
                        return path.splitext(regulator_dir[name]['file'])[0]

    return None
