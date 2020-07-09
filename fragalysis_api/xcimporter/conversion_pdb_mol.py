#!/usr/bin/env python

from rdkit import Chem
from rdkit.Chem import AllChem
import json
import os
import shutil
import warnings
import csv


class Ligand:
    def __init__(self, target_name, infile, RESULTS_DIRECTORY):
        self.infile = infile
        self.target_name = target_name
        self.mol_lst = []
        self.mol_dict = {"directory": [], "mol": [], "file_base": []}
        self.RESULTS_DIRECTORY = RESULTS_DIRECTORY
        self.non_ligs = json.load(
            open(os.path.join(os.path.dirname(__file__), "non_ligs.json"), "r")
        )
        self.pdbfile = open(os.path.abspath(infile)).readlines()
        self.hetatms = []
        self.conects = []
        self.final_hets = []
        self.wanted_ligs = []
        self.new_lig_name = "NONAME"

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
        all_ligands = []  # all ligands go in here, including solvents and ions
        for line in self.pdbfile:
            if line.startswith("HETATM"):
                all_ligands.append(line)

        for lig in all_ligands:
            if (
                lig.split()[3][-3:] not in self.non_ligs
            ):  # this takes out the solvents and ions a.k.a non-ligands
                self.wanted_ligs.append(lig[16:20].strip() + lig[20:26])
                # print(lig[16:20].strip() + lig[20:26])

        self.wanted_ligs = list(set(self.wanted_ligs))
        # print(self.wanted_ligs)

        return self.wanted_ligs

    def create_pdb_for_ligand(self, ligand, count, monomerize):
        """
        A pdb file is produced for an individual ligand, containing atomic and connection information

        params: vari pdb conversion, ligand definition, list of ligand heteroatoms and information, connection information
        returns: .pdb file for ligand
        """

        # making ligand name which identifies type of ligand and location, so is specific to the particular ligand
        ligand_name = "".join(ligand.split())

        # out directory and filename for lig pdb
        if not self.target_name in os.path.abspath(self.infile):
            if not monomerize:
                file_base = str(
                    self.target_name
                    + "-"
                    + os.path.abspath(self.infile)
                    .split("/")[-1]
                    .replace(".pdb", "")
                    .replace("_bound", "")
                    + "_"
                    + str(count)
                )
            if monomerize:
                file_base = str(self.target_name
                                + "-"
                                + os.path.abspath(self.infile)
                                .split("/")[-1]
                                .replace(".pdb", "")
                                .replace("_bound", "")
                                )
                chain = file_base.split("_")[-1]
                file_base = file_base[:-2] + "_" + str(count) + chain

        else:
            if not monomerize:
                file_base = str(
                                os.path.abspath(self.infile)
                                .split("/")[-1]
                                .replace(".pdb", "")
                                .replace("_bound", "")
                                + "_"
                                + str(count)
                            )
            if monomerize:
                file_base = str(
                                os.path.abspath(self.infile)
                                .split("/")[-1]
                                .replace(".pdb", "")
                                .replace("_bound", "")
                            )
                chain = file_base.split("_")[-1]
                file_base = file_base[:-2] + "_" + str(count) + chain

        lig_out_dir = os.path.join(self.RESULTS_DIRECTORY, file_base)

        individual_ligand = []
        individual_ligand_conect = []
        # adding atom information for each specific ligand to a list
        for atom in self.final_hets:
            if str(atom[16:20].strip() + atom[20:26]) == str(ligand):
                individual_ligand.append(atom)

        con_num = 0
        for atom in individual_ligand:
            atom_number = atom.split()[1]
            for conection in self.conects:
                if (
                    atom_number in conection
                    and conection not in individual_ligand_conect
                ):
                    individual_ligand_conect.append(conection)
                    con_num += 1

        # checking that the number of conect files and number of atoms are almost the same
        # (taking into account ligands that are covalently bound to the protein

        # assert 0 <= con_num - len(individual_ligand) <= 1

        # making into one list that is compatible with conversion to mol object
        ligand_het_con = individual_ligand + individual_ligand_conect

        # make a pdb file for the ligand molecule

        if not os.path.isdir(lig_out_dir):
            os.makedirs(lig_out_dir)

        ligands_connections = open(
            os.path.join(lig_out_dir, (file_base + ".pdb")), "w+"
        )
        for line in ligand_het_con:
            ligands_connections.write(str(line))
        ligands_connections.close()

        # making pdb file into mol object

        pdb_block = open(os.path.join(lig_out_dir, (file_base + ".pdb")), 'r').read()

        new_pdb_block = ''

        for lig in pdb_block.split('\n'):
            if 'ATM' in lig:
                pos = 16
                s = lig[:pos] + ' ' + lig[pos + 1:]
                new_pdb_block += s
            else:
                new_pdb_block += lig

            new_pdb_block += '\n'

        m = Chem.rdmolfiles.MolFromPDBBlock(new_pdb_block)

        if not m:
            print(f'WARNING: {file_base} did not produce a mol object from its pdb lig file!')

        try:
            Chem.AddHs(m)

            self.mol_lst.append(m)
            self.mol_dict["directory"].append(lig_out_dir)
            self.mol_dict["mol"].append(m)
            self.mol_dict["file_base"].append(file_base)

        except:
            print(file_base, 'is unable to produce a ligand file')

    def create_mol_file(self, directory, file_base, mol_obj, smiles_file=None):
        """
        a .mol file is produced for an individual ligand

        params: ligand definition, pdb file, pdb conversion
        returns: .mol file for the ligand
        """

        out_file = os.path.join(directory, str(file_base + ".mol"))

        if not mol_obj:
            print(f'WARNING: mol object is empty: {file_base}')

        if smiles_file:
            try:
                smiles = open(smiles_file, 'r').readlines()[0].rstrip()
                template = AllChem.MolFromSmiles(smiles)
                new_mol = AllChem.AssignBondOrdersFromTemplate(template, mol_obj)

                return Chem.rdmolfiles.MolToMolFile(new_mol, out_file)
            except Exception as e:
                print(e)
                print('failed to fit template ' + smiles_file)
                print(f'template smiles: {smiles}')
                return Chem.rdmolfiles.MolToMolFile(mol_obj, out_file)

        else:
            print(f'Warning: No smiles file: {file_base}' )

        # creating mol file
        return Chem.rdmolfiles.MolToMolFile(mol_obj, out_file)

    def create_sd_file(self, mol_obj, writer):
        """
        a molecular object defined in the pdb file is used to produce a .sdf file

        params: pdb file for the molecule, SDWriter from rdkit
        returns: .sdf file with all input molecules from each time the function is called
        """
        # creating sd file with all mol files
        return writer.write(mol_obj)

    def create_metadata_file(self, directory, file_base, mol_obj, smiles_file=None):
        """
        Metadata .csv file prepared for each ligand
        params: file_base and smiles
        returns: .mol file for the ligand
        """

        meta_out_file = os.path.join(directory, str(file_base + "_meta.csv"))
        smiles_out_file = os.path.join(directory, str(file_base + "_smiles.txt"))

        if smiles_file:
            try:
                smiles = open(smiles_file, 'r').readlines()[0].rstrip()
                # write to .txt file
                smiles_txt = open(smiles_out_file, "w+")
                smiles_txt.write(smiles)
                smiles_txt.close()

            except Exception as e:
                print(e)
                print('failed to open smiles file ' + smiles_file)

        else:
            print(f'Warning: No smiles file: {file_base}')
            # Try use mol_obj if no smiles file
            try:
                smiles = Chem.MolToSmiles(mol_obj)
            except Exception as e:
                print(e)
                print('failed to convert mol obj to smiles' + smiles_file)

                smiles = "NA"

        meta_data_dict = {'Blank':'',
                          'fragalysis_name': file_base,
                          'crystal_name': file_base.rsplit('_', 1)[0],
                          'smiles': smiles,
                          'new_smiles':'',
                          'alternate_name':'',
                          'site_name':'',
                          'pdb_entry':''}

        # Write dict to csv
        meta_data_file = open(meta_out_file, 'w+')
        w = csv.DictWriter(meta_data_file, meta_data_dict.keys())
        w.writerow(meta_data_dict)
        meta_data_file.close()

class pdb_apo:
    def __init__(self, infile, target_name, RESULTS_DIRECTORY, filebase):
        self.target_name = target_name
        self.pdbfile = open(infile).readlines()
        self.RESULTS_DIRECTORY = RESULTS_DIRECTORY
        self.filebase = filebase
        self.non_ligs = json.load(
            open(os.path.join(os.path.dirname(__file__), "non_ligs.json"), "r")
        )
        self.apo_file = None

    def make_apo_file(self):
        """
        Keeps anything other than unique ligands

        :param: pdb file
        :returns: created XXX_apo.pdb file
        """
        lines = ""
        for line in self.pdbfile:
            if (
                line.startswith("HETATM")
                and line.split()[3] not in self.non_ligs
                or line.startswith("CONECT")
            ):
                continue
            else:
                lines += line

        apo_file = open(
            os.path.join(self.RESULTS_DIRECTORY, str(self.filebase + "_apo.pdb")), "w+"
        )
        apo_file.write(str(lines))
        apo_file.close()
        self.apo_file = os.path.join(
            self.RESULTS_DIRECTORY, str(self.filebase + "_apo.pdb")
        )

    def make_apo_desol_files(self):
        """
        Creates two files:
        _apo-desolv - as apo, but without solvent, ions and buffers;
        _apo-solv - just the ions, solvent and buffers

        :returns: Created files
        """
        prot_file = open(
            os.path.join(
                self.RESULTS_DIRECTORY, str(self.filebase + "_apo-desolv.pdb")
            ),
            "w+",
        )
        solv_file = open(
            os.path.join(self.RESULTS_DIRECTORY, str(self.filebase + "_apo-solv.pdb")),
            "w+",
        )
        if not self.apo_file:
            return Warning(
                "Apo file has not been created. Use pdb_apo().make_apo_file()"
            )
        else:
            for line in open(self.apo_file).readlines():
                if line.startswith("HETATM"):
                    solv_file.write(line)
                else:
                    prot_file.write(line)
        solv_file.close()
        prot_file.close()


def set_up(target_name, infile, out_dir, monomerize, smiles_file=None):

    """

    :param pdbcode: pdb code that has already been uploaded into directory of user ID
    :param USER_ID: User ID and timestamp that has been given to user when they upload their files
    :return: for each ligand: pdb, mol files. For each pdb file: sdf and apo.pdb files.
    """

    RESULTS_DIRECTORY = os.path.join(out_dir, target_name)

    if not os.path.isdir(RESULTS_DIRECTORY):
        os.makedirs(RESULTS_DIRECTORY)

    print(RESULTS_DIRECTORY)

    new = Ligand(
        target_name, infile, RESULTS_DIRECTORY
    )  # takes in pdb file and returns specific ligand files
    new.hets_and_cons()  # takes only hetatm and conect file lines from pdb file
    new.remove_nonligands()  # removes ions and solvents from list of ligands
    new.find_ligand_names_new()  # finds the specific name and locations of desired ligands
    for i in range(len(new.wanted_ligs)):
        new.create_pdb_for_ligand(
            new.wanted_ligs[i], count=i, monomerize=monomerize
        )  # creates pdb file and mol object for specific ligand

    for i in range(len(new.mol_dict["directory"])):

        if not new.mol_dict["mol"][i]:
            warnings.warn(
                str(
                    "RDkit mol object for "
                    + new.mol_dict["file_base"][i]
                    + " is None, please check the input. Will not write any files"
                )
            )
            continue

        shutil.copy(
            infile,
            os.path.join(
                new.mol_dict["directory"][i],
                str(new.mol_dict["file_base"][i] + "_bound.pdb"),
            ),
        )

        new.create_mol_file(
            directory=new.mol_dict["directory"][i],
            file_base=new.mol_dict["file_base"][i],
            mol_obj=new.mol_dict["mol"][i],
            smiles_file=smiles_file,
        )  # creates mol file for each ligand

        writer = Chem.rdmolfiles.SDWriter(
            os.path.join(
                new.mol_dict["directory"][i],
                str(new.mol_dict["file_base"][i] + ".sdf"),
            )
        )

        new.create_sd_file(
            new.mol_lst[i], writer
        )  # creates sd file containing all mol files
        writer.close()  # this is important to make sure the file overwrites each time

        new.create_metadata_file(
            mol_obj=new.mol_dict["mol"][i],
            directory=new.mol_dict["directory"][i],
            file_base=new.mol_dict["file_base"][i],
            smiles_file=smiles_file,
        )  # create metadata csv file for each ligand

        new_apo = pdb_apo(
            infile,
            target_name,
            new.mol_dict["directory"][i],
            new.mol_dict["file_base"][i],
        )
        new_apo.make_apo_file()  # creates pdb file that doesn't contain any ligand information
        new_apo.make_apo_desol_files()  # makes apo file without solvent, ions and buffers, and file with just those

    return new
