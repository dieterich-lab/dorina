# vim: set fileencoding=utf-8 :

import logging
import os
from os import path
from dorina import config

def get_genomes():
    """Get all available genomes"""
    options = config.get_config()
    genomes = {}

    root = path.join(options.data.path, 'genomes')

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

                for gff_file in os.listdir(assembly_path):
                    gff_path = path.join(assembly_path, gff_file)
                    if not path.isfile(gff_path):
                        logging.debug("skipping non-file %r" % gff_path)
                        continue

                    root, ext = path.splitext(gff_file)
                    if ext in ('.gff', '.bed'):
                        assembly_dict[root] = True

    return genomes
