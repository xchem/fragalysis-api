#!/usr/bin/env python

from rdkit import Chem
import json
import os

class Ligand:
    def __init__(self, pdbcode_):
        self.pdbcode = pdbcode_
        self.mol_lst = []
        self.RESULTS_DIRECTORY = "../results/" + str(self.pdbcode)
        self.non_ligs = json.load(open("../non_ligs", "r"))
        self.pdbfile = open("../anna/" + str(self.pdbcode) + ".pdb").readlines()
        self.hetatms = []
        self.conects = []
        self.final_hets = []
        self.wanted_ligs = []
        self.new_lig_name = 'NONAME'

    def make_directory(self):
        try:
            os.makedirs(self.RESULTS_DIRECTORY)
        except FileExistsError:
            for file in os.listdir(self.RESULTS_DIRECTORY):
                file_path = os.path.join(self.RESULTS_DIRECTORY, file)
                os.unlink(file_path)

    def hets_and_cons(self):
        """
        Heteroatoms and connect files are pulled out from full pdb file

        param: pdb file in .readlines() format
        returns: lits of hetatomic information and connection information
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
        retunns: list of heteroatoms that are not contained in the non_ligs list
        """

        for line in self.hetatms:
            if str(line.split()[3]) not in self.non_ligs:
                self.final_hets.append(line)
        return self.final_hets

    def find_ligand_names_new(self):
        """
        Finds list of ligands contained in the structure, including
        """
        all_ligands = []
        for line in self.pdbfile:
            if line.startswith("HET "):
                all_ligands.append(line)

        for lig in all_ligands:
            if lig.split()[1] not in self.non_ligs:
                self.wanted_ligs.append(lig.split()[1:])
        return self.wanted_ligs

    def determine_ligands(self):
        self.unique = []

        for i in range(len(self.wanted_ligs)):
            if self.wanted_ligs[i][0] not in self.unique:
                self.unique.append(self.wanted_ligs[i][0])
        # if len(unique) > 1:
        #     print("More than one possible ligand has been found. These are: ")
        #     for i in unique:
        #         print(i)
        #     new_lig_name = input(
        #         "Which ligand would you like to extract?\nIf you do not have this information, leave blank.\n")
        # new_lig_info = []
        # if len(new_lig_name) == 3:
        #     for i in ligands:
        #         if new_lig_name in i:
        #             new_lig_info.append(i)
        #     return new_lig_info
        # else:
        return self.unique

    def create_pdb_for_ligand(self, ligand):
        """
        A pdb file is produced for an individual ligand, containing atomic and connection information

        params: vari pdb conversion, ligand definition, list of ligand heteroatoms and information, connection information
        returns: .pdb file for ligand
        """
        if len(ligand) == 4:
            ligand_name = str(ligand[0] + '_' + str(ligand[1]) + '_' + str(ligand[2]))
        elif len(ligand) == 3:
            ligand_name = str(ligand[0] + '_' + str(ligand[1]))
        else:
            ligand_name = str(ligand)

        individual_ligand = []
        individual_ligand_conect = []

        if len(ligand) == 3:
            for atom in self.final_hets:
                atom_new = atom.replace('HETATM', '')
                if len(atom_new.split()[1]) <= 3 and ligand[0] in atom_new.split()[2] and ligand[1] == atom_new.split()[3]:
                    individual_ligand.append(atom)
                elif len(atom_new.split()[1]) > 3 and ligand[0] in atom_new.split()[1] and ligand[1] == \
                        atom_new.split()[2]:
                    individual_ligand.append(atom)
        if len(ligand) == 4:
            for atom in self.final_hets:
                atom_new = atom.replace('HETATM', '')
                if len(atom_new.split()[1]) <= 3 and ligand[0] in atom_new.split()[2] and ligand[1] == atom_new.split()[3] and ligand[2] == atom_new.split()[4]:
                    individual_ligand.append(atom)
                elif len(atom_new.split()[1]) > 3 and ligand[0] in atom_new.split()[1] and ligand[1] == \
                        atom_new.split()[2] and ligand[2] == atom_new.split()[3]:
                    individual_ligand.append(atom)

        assert (len(individual_ligand) == int(ligand[-1]))

        con_num = 0
        for atom in individual_ligand:
            atom_new = atom.replace('HETATM', '')
            atom_number = atom_new.split()[0]
            for conection in self.conects:
                if atom_number in conection and conection not in individual_ligand_conect:
                    individual_ligand_conect.append(conection)
                    con_num += 1

        assert (0 <= con_num - len(individual_ligand) <= 1)

        ligand_het_con = individual_ligand + individual_ligand_conect

        ligands_connections = open(os.path.join(self.RESULTS_DIRECTORY, str(self.pdbcode) + "_" + str(ligand_name) + ".pdb"),"w+")
        for line in ligand_het_con:
            ligands_connections.write(str(line))
        ligands_connections.close()

        m = Chem.rdmolfiles.MolFromPDBFile(
            "../results/" + str(self.pdbcode) + "/" + str(self.pdbcode) + "_" + str(ligand_name) + ".pdb")
        self.mol_lst.append(m)

    def create_mol_file(self, ligand, mol_obj):
        """
        a .mol file is produced for an individual ligand

        params: ligand definition, pdb file, pdb conversion
        returns: .mol file for the ligand
        """
        if len(ligand) == 4:
            ligand_name = str(ligand[0] + '_' + str(ligand[1]) + '_' + str(ligand[2]))
        elif len(ligand) == 3:
            ligand_name = str(ligand[0] + '_' + str(ligand[1]))
        else:
            ligand_name = str(ligand)

        return Chem.rdmolfiles.MolToMolFile(mol_obj, os.path.join(self.RESULTS_DIRECTORY, str(self.pdbcode) + "_" + str(ligand_name) + "_mol.mol"))

    def create_sd_file(self, mol_obj, writer):
        """
        a molecular object defined in the pdb file is used to produce a .sdf file

        params: pdb file for the molecule, SDWriter from rdkit
        returns: .sdf file with all input molecules from each time the function is called
        """
        return writer.write(mol_obj)


def set_up(pdbcode):
    new = Ligand(pdbcode)
    new.make_directory()
    new.hets_and_cons()
    new.remove_nonligands()
    new.find_ligand_names_new()
    if len(new.wanted_ligs) > 1:
        new.determine_ligands()
    for i in range(len(new.wanted_ligs)):
        new.create_pdb_for_ligand(new.wanted_ligs[i])

    writer = Chem.rdmolfiles.SDWriter(os.path.join(new.RESULTS_DIRECTORY, str(new.pdbcode) + "_out.sdf"))
    for i in range(len(new.mol_lst)):
        new.create_mol_file(new.wanted_ligs[i], new.mol_lst[i])
        new.create_sd_file(new.mol_lst[i], writer)
    writer.close()  # this is important to make sure the file overwrites

    return new

set_up('2bui')