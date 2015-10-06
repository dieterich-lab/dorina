import os
import re

class Genome:
    @staticmethod
    def parse_func(root, assembly_dict):
        """Parse function used to initialise all genomes from the data directory.
        """

        for gff_file in os.listdir(root):
            gff_path = os.path.join(root, gff_file)
            if not os.path.isfile(gff_path):
                continue

            basename, ext = os.path.splitext(gff_file)
            if ext in ('.gff', '.bed'):
                assembly_dict[basename] = True

        return assembly_dict

    # TODO: turn this into a regular method on genomes
    @staticmethod
    def get_genes(name):
        """Get a list of genes from genome <name>"""
        genes = []

        # TODO: only do this once when initialising the genome instance
        gene_name = re.compile(r'.*ID=(.*?)($|;\w+)')

        genome_dir = Genome.path_by_name(name)
        if genome_dir is None:
            return genes

        genome = os.path.join(genome_dir, 'all.gff')
        if not os.path.exists(genome):
            return genes

        for line in open(genome, 'r'):
            match = gene_name.match(line)
            if match is not None:
                genes.append(match.group(1))

        return genes
