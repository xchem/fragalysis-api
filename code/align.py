import glob
#import Bio.PDB
from pymol import cmd
import os

class Align:
    def __init__(self, directory):
        self.directory = directory
        self.ref = None
        self.load_obj = None
        self.aligned_obj = None

    @property
    def pdb_in_list(self):
        return sorted(glob.glob(os.path.join(self.directory, "*.pdb")))


    def load(self):
        for i in self.pdb_in_list:
            cmd.load()
        self.load_obj = cmd


    def pick_ref(self):
        pass

    def align(self):
        self.aligned_obj = cmd.align(self.load_obj, self.ref)

    # method of splitting the file projeduce from alignment.

    def save_pymol_objs(self):
        """ Method to save a list of objects as individual pdb files.
        """
        #cmd.multifilesave(/tmp/{self.align_obs}.pdb)
        #cmd.multifilesave()
    def split_merged(self):
        """ Method to split a pdb file containing multiple structures
            into individual objects, before calling another method to
            save each object as individual pdb files.
        """
        cmd.split_states(self.merged_align , prefix = self.directory)
        cmd.delete(self.merged_align)
        self.align_obs = cmd.get_object_list('(all)')
        Align.save_pymol_objs()