from pymol import cmd


class Align:
    def __init__(self, directory="", pdb_in_list=[]):
        self.reference = []
        self.directory = directory
        self.pdb_in_list = self.gen_list()
        self.align = []
        self.pdb_out_list = []
        self.align_obs = []
        self.merged_align = []

    def align(self):
        for i in self.pdb_in_list:
            cmd.load(i)
            self.merged_align = cmd.align(self.pdb_in_list(i), self.reference)
        return self.merged_align

