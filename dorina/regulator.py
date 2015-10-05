import os
from pybedtools import BedTool

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

    @staticmethod
    def merge(regulators):
        """Merge a list of regulators using BedTool.cat"""
        if len(regulators) > 1:
            return BedTool.cat(*regulators, postmerge=False)
        else:
            return regulators[0]

    @staticmethod
    def from_name(directory, name_or_path, assembly=None):
        if os.sep in name_or_path:
            return Regulator("custom", name_or_path, True)

        if not assembly:
            raise ValueError("Must provide assembly")

        filename = None
        for species, species_dir in directory.items():
            for _assembly, assembly_dir in species_dir.items():
                if assembly == _assembly and name_or_path in assembly_dir:
                    basename = os.path.splitext(assembly_dir[name_or_path]['file'])[0]
                    filename = basename + ".bed"

        if not filename or not os.path.isfile(filename):
            raise ValueError("Could not find regulator: %s" % name_or_path)

        return Regulator(name_or_path, filename, False)

    @staticmethod
    def from_names(directory, names, assembly):
        if names:
            return map(lambda x: Regulator.from_name(directory, x, assembly).bed,
                       names)
        else:
            return []
