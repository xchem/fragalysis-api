import glob
import Bio.PDB as bp
import pymol
import pandas as pd
import os
import warnings
import Bio.PDB.PDBExceptions as bpp
import json
import shutil
warnings.simplefilter('ignore', bpp.PDBConstructionWarning)


class Align:

    def __init__(self, directory, pdb_ref=''):

        self.directory = directory
        self._get_ref = pdb_ref

    @property
    def _get_files(self):
        """
        Extracts a list of paths for all PDBs within the given directory.
        :return: list of .pdb file names in directory
        """
        all_files = set(glob.glob(os.path.join(self.directory, '*')))
        txt_files = set(glob.glob(os.path.join(self.directory, '*.txt')))
        return list(all_files-txt_files)

    def _load_objs(self):
        """
        Loads each pdb into the PyMol instance/object.
        :type: object
        :return: PyMol object with .pdb protein structure file loaded
        """
        # Looping through each pdb file in the directory and loading them into the cmd
        for num, file in enumerate(self._get_files):
            pymol.cmd.load(file, os.path.splitext(os.path.basename(file))[0])

        # deal with files that have no conect records
        pymol.cmd.set('pdb_conect_all', 'on')

        return pymol.cmd

    @property
    def _get_ref(self):
        """
        Determines the best reference structure for alignments if not user provided.
        Chosen based on longest length with lowest resolution.
        :return: str of .pdb filename of reference:
        """
        return self.__pdb_ref

    @_get_ref.setter
    def _get_ref(self, pdb_ref):
        """
        Determines the best reference structure for alignments if not user provided.
        Chosen based on longest length with lowest resolution.
        :param pdb_ref: pdb to use as reference pdb for alignments
        :return PyMol instance with reference object assigned as reference property:
        """
        if pdb_ref != '':
            try:
                assert(os.path.isfile(os.path.join(self.directory, pdb_ref+".pdb")) == 1)
                self.__pdb_ref = pdb_ref
            except AssertionError:
                print('pdb desired as reference does not exists. Default pdb chosen.')
                self.__pdb_ref = self.__best_length_and_resolution([i for i in self._get_files if 'pdb' in i])
        else:
            self.__pdb_ref = self.__best_length_and_resolution([i for i in self._get_files if 'pdb' in i])

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
        Determine resolution, sequence and length of .pdb file.
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

    def __get_pdb_file(self, pdb_name):
        """
        Returns the pdb file path for a given pdb name.
        :param pdb_name: name of pdb
        :return: file path to the pdb with the given name
        """

        for file in self._get_files:
            if pdb_name == os.path.splitext(os.path.basename(file))[0]:
                return file

    def __get_header(self, pdb_file):
        """
        Identifies the section of a PDB which contains the headers ATOM/HETATM
        :param pdb_file: The pdb to acquire the header locations of
        :return: front locations of the ATOM/HETATM headers in the given pdb, end locations of the ATOM/HETATM headers in the given pdb       
        """
        with open(pdb_file) as handle:
            switch = 0
            header_front, header_end = [], []

            for line in handle:

                if line.startswith('ATOM'): switch = 1

                if line.startswith('HETATM'): switch = 2

                if switch == 0: header_front.append(line)

                if (switch == 2) and not line.startswith('HETATM'): header_end.append(line)

        return header_front, header_end

    def _save_align(self, name, pdb_base, out_dir):
        """
        Saves modified pdb as pdb file again. It also ensures the pdb header is kept in the new file.
        :param name: name of pdb
        :param pdb_base: coordinates of atom in the pdb format
        :param out_dir: directory to save new pdb file in
        :return: a saved pdb file
        """
        if not os.path.exists(out_dir):  # Creating output directory if it doesn't already exist
            os.makedirs(out_dir)

        pdb_file = self.__get_pdb_file(name)
        header_front = self.__get_header(pdb_file)[0]

        with open(os.path.join(out_dir, f'{name}_bound.pdb'), 'w') as handle:

            new_pdb = header_front + [pdb_base]

            for line in new_pdb:
                handle.write(line)

    def align(self, out_dir):
        """
        Aligns all pdbs in the pymol object to the pdb_ref.
        :param out_dir: directory to save aligned pdbs in
        :return: saves the pdbs
        """

        # Silently open PyMOL
        pymol.pymol_argv = ['pymol', '-qc']
        pymol.finish_launching()
        pymol_cmd = self._load_objs()

        for num, name in enumerate(pymol_cmd.get_names()):  # Saves the aligned pdb files from the cmd as pdb files

            if not name == self._get_ref:
                pymol_cmd.align(name, self._get_ref)
                self._save_align(name, pymol_cmd.get_pdbstr(selection=name), out_dir)

            elif name == self._get_ref:
                self._save_align(name, pymol_cmd.get_pdbstr(selection=name), out_dir)


class Monomerize:

    def __init__(self, directory, outdir):
        self.directory = directory
        self.outdir = outdir
        self.non_ligs = json.load(
            open(os.path.join(os.path.dirname(__file__), "non_ligs.json"), "r")
        )

    def get_filelist(self):
        return glob.glob(os.path.join(self.directory, '*.pdb'))

    def find_ligs(self, pdb_lines):
        """
        Finds list of ligands contained in the structure, including
        """
        all_ligands = []  # all ligands go in here, including solvents and ions
        wanted_ligs = []
        for line in pdb_lines:
            if line.startswith("HETATM"):
                all_ligands.append(line)

        for lig in all_ligands:
            if (
                    lig.split()[3][-3:] not in self.non_ligs
            ):  # this takes out the solvents and ions a.k.a non-ligands
                wanted_ligs.append(lig[16:20] + lig[20:26])
                # print(lig[16:20].strip() + lig[20:26])

        wanted_ligs = list(set(wanted_ligs))

        return wanted_ligs

    def save_chain(self, lig, f):
        lig_chain = lig[5]
        name = os.path.splitext(os.path.basename(f))[0] + '_' + str(lig_chain)
        filename = os.path.join(self.outdir, f'{name}_mono.pdb')
        pymol.cmd.load(f, name)
        pymol.cmd.save(filename, f'chain {lig_chain}')
        pymol.cmd.reinitialize()

        return filename

    def process_ligs(self, filename):
        test_block = open(filename, 'r').readlines()
        ligs = self.find_ligs(test_block)
        print(ligs)
        outnames = []
        for lig in ligs:
            o = self.save_chain(lig, filename)
            outnames.append(o)
            if os.path.isfile(filename.replace('.pdb', '_smiles.txt')):
                shutil.copy(filename.replace('.pdb', '_smiles.txt'), o.replace('_mono.pdb', '_smiles.txt'))
        return outnames

    def write_bound(self, inname, outname):
        with open(inname, 'r') as handle:
            switch = 0
            header_front, header_end = [], []

            for line in handle:

                if line.startswith('ATOM'): switch = 1

                if line.startswith('HETATM'): switch = 2

                if switch == 0 and not line.startswith('REMARK 350'): header_front.append(line)

                if (switch == 2) and not line.startswith('HETATM'): header_end.append(line)

        for o in outname:
            newfile_contents = open(o, 'r').readlines()

            with open(o.replace('_mono.pdb', '.pdb'), 'w') as handle:
                remark = ['REMARK warning: chains may be ommitted for alignment\n']
                new_pdb = ''.join(remark + header_front + newfile_contents)
                print(new_pdb)
                handle.write(new_pdb)

    def monomerize_all(self):
        for f in self.get_filelist():
            print(f)
            outnames = self.process_ligs(f)
            print(outnames)
            self.write_bound(f, outnames)
            for o in outnames:
                if os.path.isfile(o):
                    os.remove(o)
