# vim: set fileencoding=utf-8 :

import logging
import os
import json

class DorinaUtils:
    @staticmethod
    def walk_assembly_tree(root, parse_func):
        """Walk a directory structure containg clade, species, assembly

        Call parse_func() for every assembly directory"""
        genomes = {}

        for species in os.listdir(root):
            species_path = os.path.join(root, species)
            if not os.path.isdir(species_path):
                #logging.debug("skipping non-directory %r" % species_path)
                continue

            description_file = os.path.join(species_path, 'description.json')
            species_dict = {}
            try:
                with open(description_file, 'r') as fh:
                    genomes[species] = json.load(fh)
                    genomes[species]['assemblies'] = species_dict
            except IOError:
                genomes[species] = species_dict

            for assembly in os.listdir(species_path):
                assembly_path = os.path.join(species_path, assembly)
                if not os.path.isdir(assembly_path):
                    #logging.debug("skipping non-directory %r" % assembly_path)
                    continue

                assembly_dict = {}
                species_dict[assembly] = assembly_dict

                parse_func(assembly_path, assembly_dict)

        return genomes
