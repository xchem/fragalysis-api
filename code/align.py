import glob
import Bio.PDB as bp
from pymol import cmd
import pandas as pd
import os


class Align:
    def __init__(self, directory, pdb_ref=''):
        self.directory = directory
        self.get_ref = pdb_ref

        # self.ref = self._pick_ref(pdb_ref)
        # self.load_obj = ""
        # self.aligned_obj = ""

    @property
    def get_files(self):
        """Extracts a list of paths for all pdbs within the given directory.

        Returns:
            list: list containing the path for each pdb in the given directory
        """
        return glob.glob(os.path.join(self.directory, "*.pdb"))

    def __load_objs(self):
        """Loads each pdb into the cmd object.

        Returns:
            object: object containing each pdb.
        """
        for i in self.get_files:
            cmd.load(i)
        return cmd

    @property
    def get_ref(self):
        """Find the pdb to use as reference for alignments. If none given as input,
        find one based on the longest length with the lowest resolution.

        Returns:
            str: pdb to use as reference
        """

        return self.__pdb_ref

    @get_ref.setter
    def get_ref(self, pdb_ref):

        if pdb_ref != '':
            # assert(ref in directory)
            self.__pdb_ref = pdb_ref
        else:
            self.__pdb_ref = self.__best_length_and_resolution(self.get_files)

    def __best_length_and_resolution(self, pdb_files):
        """Find the longest pdb structure with the lowest resolution.

        Args:
            list of files (list): filenames all pdbs

        Returns:
            str: str containing the name of the pdb chosen as reference.

        """
        a_df = pd.DataFrame(pdb_files, columns=['file'])
        a_df[['resolution', 'p_len', 'pdb']] = a_df.file.apply(self.__get_length_and_resolution)

        return a_df.sort_values(by=['p_len', 'resolution'], ascending=[False, True]).pdb[0]

    def __get_length_and_resolution(self, file):
        """Find the resolution, sequence length and pdb name for each pdb

        Args:
            file (str): filename of a single pdb

        Returns:
            pandas series: series with resolution, sequence length and pdb name
        """
        parser = bp.PDBParser()
        ppb = bp.PPBuilder()
        structure = parser.get_structure(os.path.splitext(os.path.basename(file))[0], file)
        pp = ppb.build_peptides(structure)[0]

        return pd.Series([structure.header['resolution'], len(pp.get_sequence()), structure.id])

    @property
    def get_alignment(self):
        print(os.path.join(self.directory, self.get_ref + ".pdb"))
        return cmd.align(self.__load_objs, os.path.join(self.directory, self.get_ref + ".pdb"))

    # method of splitting the file projeduce from alignment.

    def save_pymol_objs(self):
        """ Method to save a list of objects as individual pdb files.

        """
        # cmd.multifilesave(/tmp/{self.align_obs}.pdb)
        # cmd.multifilesave()

    def split_merged(self):
        """ Method to split a pdb file containing multiple structures
            into individual objects, before calling another method to
            save each object as individual pdb files.
        """
        cmd.split_states(self.merged_align, prefix=self.directory)
        cmd.delete(self.merged_align)
        self.align_obs = cmd.get_object_list('(all)')
        Align.save_pymol_objs()


struc = Align("../data/ATAD2", pdb_ref='')
print(struc.get_files)
print(struc.get_ref)
# print(struc.pdb_objs)
print(struc.get_alignment)
# print(struc.pdb_multi_saved)
# print(struc.split_merged())