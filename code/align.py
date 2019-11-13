import glob
import Bio.PDB as bp
from pymol import cmd
import pandas as pd
import os


class Align:

    def __init__(self, directory, pdb_ref=''):

        self.directory = directory
        self._get_ref = pdb_ref

    @property
    def _get_files(self):
        """
        Extracts a list of paths for all pdbs within the given directory.
        :param self:
        :type self:
        :return list of .pdb filenames in directory:
        """

        return glob.glob(os.path.join(self.directory, "*.pdb"))

    def _load_objs(self):
        """
        Loads each pdb into the PyMol instance/object.
        :param self:
        :type object:
        :return PyMol object with .pdb protein structure file loaded:
        """
        for num, file in enumerate(self._get_files):
            cmd.load(file, os.path.splitext(os.path.basename(file))[0])

        return cmd

    @property
    def _get_ref(self):
        """
        Determines the best reference structure for alignments if not user provided. Chosen based on longest length with lowest resolution.
        :param self:
        :type: object:
        :return: str of .pdb filename of reference:
        """
        return self.__pdb_ref

    @_get_ref.setter
    def _get_ref(self, pdb_ref):
        """
        :param pdb_ref:
        :type str:
        :return PyMol instance with reference object assigned as reference property:
        """
        if pdb_ref != '':
            # assert(ref in directory)
            self.__pdb_ref = pdb_ref
        else:
            self.__pdb_ref = self.__best_length_and_resolution(self._get_files)

    def __best_length_and_resolution(self, pdb_files):
        """
        Find the longest pdb structure with the lowest resolution from all imported files.
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
        pp = ppb.build_peptides(structure)[0]

        return pd.Series([structure.header['resolution'], len(pp.get_sequence()), structure.id])

    def _save_align(self):
        """
        Saves aligned structures as .pdb files
        :param self:
        :return aligned structures saved as .pdb files:
        """
        pymol_cmd = self._load_objs()

        if not os.path.exists('../data/aligned'):
            os.makedirs('../data/aligned')

        for num, name in enumerate(pymol_cmd.get_names()):

            pymol_cmd.align(name, self._get_ref)
            pymol_cmd.save(f'../data/aligned/{name}_aligned.pdb', name)
    
    def align(self):
        """
        A single method that calls the methods in sequence required to align
        the pdb files of the structure. 
        """
        self._get_files
        self._save_align()
        