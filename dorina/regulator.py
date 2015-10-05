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
