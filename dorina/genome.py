import os

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
