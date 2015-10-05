import os

class Genome:
    @staticmethod
    def path_by_name(directory, name):
        """Take a genome name and return the path to the genome directory"""
        for species, species_dir in directory.items():
            if name in species_dir['assemblies']:
                # TODO: the path should be recorded somewhere
                filename = os.path.join(self.datadir, 'genomes', species, name)

        if not filename or not os.path.isfile(filename):
            raise ValueError("Could not find genome: %s" % name)

        return filename

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
