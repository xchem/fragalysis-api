import glob
import Bio.PDB as bp
import pymol
import pandas as pd
import os
import warnings
import Bio.PDB.PDBExceptions as bpp
warnings.simplefilter('ignore', bpp.PDBConstructionWarning)


class Align:

    def __init__(self, directory, pdb_ref=''):

        self.directory = directory
        self._get_ref = pdb_ref

    @property
    def _get_files(self):
        """
        Extracts a list of paths for all PDBs within the given directory.
        :param self:
        :type self:
        :return list of .pdb file names in directory:
        """

        return glob.glob(os.path.join(self.directory, "*.pdb"))

    def _load_objs(self):
        """
        Loads each pdb into the PyMol instance/object.
        :param self:
        :type object:
        :return PyMol object with .pdb protein structure file loaded:
        """
        # Looping through each pdb file in the directory and loading them into the cmd
        for num, file in enumerate(self._get_files):
            pymol.cmd.load(file, os.path.splitext(os.path.basename(file))[0])

        return pymol.cmd

    @property
    def _get_ref(self):
        """
        Determines the best reference structure for alignments if not user provided.
        Chosen based on longest length with lowest resolution.
        :param self:
        :type: object:
        :return: str of .pdb filename of reference:
        """
        return self.__pdb_ref

    @_get_ref.setter
    def _get_ref(self, pdb_ref):
        """
        Determines the best reference structure for alignments if not user provided.
        Chosen based on longest length with lowest resolution.
        :param pdb_ref:
        :type str:
        :return PyMol instance with reference object assigned as reference property:
        """
        if pdb_ref != '':
            try:
                assert(os.path.isfile(os.path.join(self.directory, pdb_ref+".pdb")) == 1)
                self.__pdb_ref = pdb_ref
            except AssertionError:
                print('pdb desired as reference does not exists. Default pdb chosen.')
                self.__pdb_ref = self.__best_length_and_resolution(self._get_files)
        else:
            self.__pdb_ref = self.__best_length_and_resolution(self._get_files)

    def __best_length_and_resolution(self, pdb_files):
        """
        Find the longest pdb structure with the lowest resolution from all imported files. 
        This will be used as the reference pdb for alignment against.
        :param pdb_files:
        :return: str of filename with best .pdb file
        """
        a_df = pd.DataFrame(pdb_files, columns=['file'])
        a_df[['resolution', 'p_len', 'pdb']] = a_df.file.apply(self.__get_length_and_resolution)

        return a_df.sort_values(by=['p_len', 'resolution'], ascending=[False, True]).pdb[0]

    def __get_length_and_resolution(self, file):
        """
        Determine resolution, sequence and length of .pdb file
        :param file:
        :return pandas series with resolution, sequence and length and .pdb filename:
        """
        parser = bp.PDBParser()
        ppb = bp.PPBuilder()
        structure = parser.get_structure(os.path.splitext(os.path.basename(file))[0], file)

        seq_len = 0

        for pp in ppb.build_peptides(structure):  # Retrieve length by looping through each chain in the protein
            seq_len += len(pp.get_sequence())

        # using a functions from PDBParser parser class to get the resolution and protein id from the pdb file 
        return pd.Series([structure.header['resolution'], seq_len, structure.id])

    def _save_align(self, path_save):
        """
        Saves aligned structures as .pdb files
        :param self:
        :return aligned structures saved as .pdb files:
        """
        # Silently open PyMOL
        pymol.pymol_argv = ['pymol', '-qc']
        pymol.finish_launching()
        pymol_cmd = self._load_objs()

        if not os.path.exists(path_save):  # Creating output directory if it doesn't already exist
            os.makedirs(path_save)

        for num, name in enumerate(pymol_cmd.get_names()):  # Saves the aligned pdb files from the cmd as pdb files
            if not name == self._get_ref:
                pymol_cmd.align(name, self._get_ref)
                pymol_cmd.save(os.path.join(path_save, f'{name}_bound.pdb'), name)

            elif name == self._get_ref:
                pymol_cmd.save(os.path.join(path_save, f'{name}_bound.pdb'), name)

    def align(self, out_dir):
        """
        A single method that calls the methods in sequence required to align
        the pdb files of the structure.
        """
        self._save_align(out_dir)
