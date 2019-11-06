import glob
import Bio.PDB as bp
from pymol import cmd
import pandas as pd
import os


class Align:
    def __init__(self, directory, pdb_ref=''):
        self.directory = directory
        self.get_ref = pdb_ref

    @property
    def get_files(self):
        """Extracts a list of paths for all pdbs within the given directory.

        Returns:
            list: list containing the path for each pdb in the given directory
        """
        return glob.glob(os.path.join(self.directory, "*.pdb"))

    def __load_objs(self):
        """Loads each pdb into the pymol instance.

        Returns:
            object: pymol cmd object containing each pdb.
        """
        for num, file in enumerate(struc.get_files):
            cmd.load(file, os.path.splitext(os.path.basename(file))[0])

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

    def save_align(self):
        """Aligns the pdbs and saves them individually.

        :return: Saved aligned pdbs
        """
        pymol_cmd = self.__load_objs()

        for num, name in enumerate(pymol_cmd.get_names()):

            pymol_cmd.align(name, self.get_ref)
            pymol_cmd.save(f'../data/aligned/mol_{num}.pdb', name)



if __name__ == "__main__":
    struc = Align("../data/ATAD2", pdb_ref='')
    print(struc.get_files)
    print(struc.get_ref)
    struc.save_align()