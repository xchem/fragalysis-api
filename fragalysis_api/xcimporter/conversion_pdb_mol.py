#!/usr/bin/env python

from rdkit import Chem
import json
import os


class Ligand:
    def __init__(self, pdbcode_, DATA_DIRECTORY_INPUT, RESULTS_DIRECTORY):
        self.pdbcode = pdbcode_
        self.mol_lst = []
        self.RESULTS_DIRECTORY = RESULTS_DIRECTORY
        self.non_ligs = json.load(open(os.path.join(os.path.dirname(__file__), "non_ligs.json"), "r"))
        self.pdbfile = open(os.path.join(DATA_DIRECTORY_INPUT, str(self.pdbcode) + ".pdb")).readlines()
        self.hetatms = []
        self.conects = []
        self.final_hets = []
        self.wanted_ligs = []
        self.new_lig_name = 'NONAME'

    def hets_and_cons(self):
        """
        Heteroatoms and connect files are pulled out from full pdb file

        param: pdb file in .readlines() format
        returns: lists of hetatomic information and connection information
        """

        for line in self.pdbfile:
            if line.startswith("HETATM"):
                self.hetatms.append(line)
            if line.startswith("CONECT"):
                self.conects.append(line)

        return self.hetatms, self.conects

    def remove_nonligands(self):
        """
        Non-ligands such as solvents and ions are removed from the list of heteroatoms,
        to ideally leave only the target ligand

        params: list of heteroatoms and their information, list of non-ligand small molecules
            that could be in crystal structure
        returns: list of heteroatoms that are not contained in the non_ligs list
        """

        for line in self.hetatms:
            ligand_name = line[17:20].strip()
            if ligand_name not in self.non_ligs:
                self.final_hets.append(line)
        return self.final_hets

    def find_ligand_names_new(self):
        """
        Finds list of ligands contained in the structure, including
        """
        all_ligands = []   # all ligands go in here, including solvents and ions
        for line in self.pdbfile:
            if line.startswith("HET "):
                all_ligands.append(line)

        for lig in all_ligands:
            if lig.split()[1] not in self.non_ligs: # this takes out the solvents and ions a.k.a non-ligands
                self.wanted_ligs.append(lig.split()[1:])
        return self.wanted_ligs

    def create_pdb_for_ligand(self, ligand):
        """
        A pdb file is produced for an individual ligand, containing atomic and connection information

        params: vari pdb conversion, ligand definition, list of ligand heteroatoms and information, connection information
        returns: .pdb file for ligand
        """

        # making ligand name which identifies type of ligand and location, so is specific to the particular ligand
        if len(ligand) == 4:
            ligand_name = str(ligand[0] + '_' + str(ligand[1]) + '_' + str(ligand[2]))
        elif len(ligand) == 3:
            ligand_name = str(ligand[0] + '_' + str(ligand[1]))
        else:
            ligand_name = str(ligand)

        individual_ligand = []
        individual_ligand_conect = []
        # adding atom information for each specific ligand to a list
        if len(ligand) == 3:
            for atom in self.final_hets:
                if atom[17:20].strip() == ligand[0] and atom[21:26].strip() == ligand[1]:
                    individual_ligand.append(atom)

        if len(ligand) == 4:
            for atom in self.final_hets:
                if atom[17:20].strip() == ligand[0] and atom[21:22].strip() == ligand[1] and atom[22:26].strip() == ligand[2]:
                    individual_ligand.append(atom)

        assert (len(individual_ligand) == int(ligand[-1]))

        # finding the conect files for the individual ligand
        con_num = 0
        for atom in individual_ligand:
            atom_new = atom.replace('HETATM', '')
            atom_number = atom_new.split()[0]
            for conection in self.conects:
                if atom_number in conection and conection not in individual_ligand_conect:
                    individual_ligand_conect.append(conection)
                    con_num += 1

        # checking that the number of conect files and number of atoms are almost the same
        # (taking into account ligands that are covalently bound to the protein
        assert (0 <= con_num - len(individual_ligand) <= 1)

        # making into one list that is compatible with conversion to mol object
        ligand_het_con = individual_ligand + individual_ligand_conect

        # make a pdb file for the ligand molecule
        ligands_connections = open(os.path.join(self.RESULTS_DIRECTORY, str(self.pdbcode) + "_" + str(ligand_name) + ".pdb"),"w+")
        for line in ligand_het_con:
            ligands_connections.write(str(line))
        ligands_connections.close()

        # making pdb file into mol object
        m = Chem.rdmolfiles.MolFromPDBFile(
            os.path.join(self.RESULTS_DIRECTORY, str(self.pdbcode) + "_" + str(ligand_name) + ".pdb"))
        self.mol_lst.append(m)

    def create_mol_file(self, ligand, mol_obj):
        """
        a .mol file is produced for an individual ligand

        params: ligand definition, pdb file, pdb conversion
        returns: .mol file for the ligand
        """

        # making ligand name
        if len(ligand) == 4:
            ligand_name = str(ligand[0] + '_' + str(ligand[1]) + '_' + str(ligand[2]))
        elif len(ligand) == 3:
            ligand_name = str(ligand[0] + '_' + str(ligand[1]))
        else:
            ligand_name = str(ligand)

        # creating mol file
        return Chem.rdmolfiles.MolToMolFile(mol_obj, os.path.join(self.RESULTS_DIRECTORY, str(self.pdbcode) + "_" + str(ligand_name) + "_mol.mol"))

    def create_sd_file(self, mol_obj, writer):
        """
        a molecular object defined in the pdb file is used to produce a .sdf file

        params: pdb file for the molecule, SDWriter from rdkit
        returns: .sdf file with all input molecules from each time the function is called
        """
        # creating sd file with all mol files
        return writer.write(mol_obj)


class pdb_apo:
    def __init__(self, pdbcode_, DATA_DIRECTORY_INPUT, RESULTS_DIRECTORY):
        self.pdbcode = pdbcode_
        self.pdbfile = open(os.path.join(DATA_DIRECTORY_INPUT, str(self.pdbcode) + ".pdb")).readlines()
        self.RESULTS_DIRECTORY = RESULTS_DIRECTORY

    def make_apo_file(self):
        """
        Finds terminus of protein chain and takes everything up to and including terminus into a pdb file

        :param: pdb file
        :returns: created XXX_apo.pdb file
        """
        apo_file_lst = []
        line_ter = 0
        for i in range(len(self.pdbfile)):
            if self.pdbfile[i].startswith('TER'): # finds terminus in PDB file
                line_ter = i+1

        for i in range(line_ter):  # takes every line up to and including terminus and appends to list
            apo_file_lst.append(self.pdbfile[i])

        apo_file = open(os.path.join(self.RESULTS_DIRECTORY, self.pdbcode + "_apo.pdb"), "w+")
        for line in apo_file_lst:  # turns list into pdb file
            apo_file.write(str(line))
        apo_file.close()


def set_up(pdbcode, USER_ID, in_dir, out_dir):

    """

    :param pdbcode: pdb code that has already been uploaded into directory of user ID
    :param USER_ID: User ID and timestamp that has been given to user when they upload their files
    :return: for each ligand: pdb, mol files. For each pdb file: sdf and apo.pdb files.
    """

    DATA_DIRECTORY_INPUT = os.path.join(in_dir, USER_ID)
    RESULTS_DIRECTORY = os.path.join(out_dir, USER_ID, "tmp")

    if not os.path.isdir(RESULTS_DIRECTORY):
        os.makedirs(RESULTS_DIRECTORY)

    new = Ligand(pdbcode, DATA_DIRECTORY_INPUT, RESULTS_DIRECTORY)  # takes in pdb file and returns specific ligand files
    new.hets_and_cons()  # takes only hetatm and conect file lines from pdb file
    new.remove_nonligands()  # removes ions and solvents from list of ligands
    new.find_ligand_names_new()  # finds the specific name and locations of desired ligands
    for i in range(len(new.wanted_ligs)):
        new.create_pdb_for_ligand(new.wanted_ligs[i])  # creates pdb file and mol object for specific ligand

    writer = Chem.rdmolfiles.SDWriter(os.path.join(new.RESULTS_DIRECTORY, str(new.pdbcode) + "_out.sdf"))
    for i in range(len(new.mol_lst)):
        new.create_mol_file(new.wanted_ligs[i], new.mol_lst[i])  # creates mol file for each ligand
        new.create_sd_file(new.mol_lst[i], writer)  # creates sd file containing all mol files
    writer.close()  # this is important to make sure the file overwrites each time

    new_apo = pdb_apo(pdbcode, DATA_DIRECTORY_INPUT, RESULTS_DIRECTORY)
    new_apo.make_apo_file() # creates pdb file that doesn't contain any ligand information

    return new

