#!/usr/bin/env python
import glob

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
import json
import os
import shutil
import warnings
import csv
import pypdb
import numpy as np
import gemmi


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

        :return: lists of hetatomic information and connection information
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
        :return: list of heteroatoms that are not contained in the non_ligs list
        """

        for line in self.hetatms:
            ligand_name = line[17:20].strip()
            if ligand_name not in self.non_ligs:
                self.final_hets.append(line)
        return self.final_hets

    def find_ligand_names_new(self):
        """
        Finds list of ligands contained in the structure, including solvents and ions
        :return: A listed of ligands that are not listed in the non_ligs.json file!
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

    def get_3d_distance(self, coord_a, coord_b):
        '''
        Get the 3D distance between two points
        :param coord_a: List of len 3 corresponding to xyz coords
        :param coord_b: List of len 3 corresponding to xyz coords
        :return: The distance.
        '''
        sum_ = (sum([(float(coord_a[i]) - float(coord_b[i])) ** 2 for i in range(3)]))
        return np.sqrt(sum_)

    def handle_covalent_mol(self, lig_res_name, non_cov_mol):
        '''
        Do some magic if we think the molecule has a covalent attachment
        :param lig_res_name: Name of the covalent ligand
        :param non_cov_mol: Previous .mol file that does not have covalent attachment in it.
        :return: A new mol file IF the lig_res is indeed covalent.
        '''
        # original pdb = self.pdbfile (already aligned)
        # lig res name = name of ligand to find link for

        covalent = False

        for line in self.pdbfile:
            if 'LINK' in line:
                zero = line[13:27]
                one = line[43:57]

                if lig_res_name in zero:
                    res = one
                    covalent = True

                if lig_res_name in one:
                    res = zero
                    covalent = True

        if covalent:
            for line in self.pdbfile:
                if 'ATOM' in line and line[13:27] == res:
                    res_x = float(line[31:39])
                    res_y = float(line[39:47])
                    res_z = float(line[47:55])
                    res_coords = [res_x, res_y, res_z]
                    print(res_coords)
                    atm = Chem.MolFromPDBBlock(line)
                    atm_trans = atm.GetAtomWithIdx(0)

            orig_pdb_block = Chem.MolToPDBBlock(non_cov_mol)

            lig_block = '\n'.join([l for l in orig_pdb_block.split('\n') if 'COMPND' not in l])
            lig_lines = [l for l in lig_block.split('\n') if 'HETATM' in l]
            j = 0
            old_dist = 100
            for line in lig_lines:
                j += 1
                #                 print(line)
                if 'HETATM' in line:
                    coords = [line[31:39].strip(), line[39:47].strip(), line[47:55].strip()]
                    dist = self.get_3d_distance(coords, res_coords)

                    if dist < old_dist:
                        ind_to_add = j
                        print(dist)
                        old_dist = dist

            i = non_cov_mol.GetNumAtoms()
            edmol = Chem.EditableMol(non_cov_mol)
            edmol.AddAtom(atm_trans)
            edmol.AddBond(ind_to_add - 1, i, Chem.BondType.SINGLE)
            new_mol = edmol.GetMol()
            conf = new_mol.GetConformer()
            conf.SetAtomPosition(i, Point3D(res_coords[0], res_coords[1], res_coords[2]))

            return new_mol

    def create_pdb_mol(self, file_base, lig_out_dir, smiles_file, handle_cov=False):
        """
        :param file_base: fragalysis crystal name
        :param lig_out_dir: output directory
        :param smiles_file: smiles file associated with pdb
        :param handle_cov: bool to indicate if output mol file should account of
                covalent attachment to model
        :return: mol object that attempts to correct bond order if PDB entry
                or mol object extracted from pdb file
        """
        pdb_block = open(os.path.join(lig_out_dir, (file_base + ".pdb")), 'r').read()
        print('mol1')
        lig_line = open(os.path.join(lig_out_dir, (file_base + ".pdb")), 'r').readline()
        res_name = lig_line[16:20].replace(' ', '')
        print('mol2')
        print(res_name)
        # Look for PDB entries in PDB bank and use residue name to get bond order
        if not smiles_file:
            print('mol3')
            try:
                print('mol4')
                print(pdb_block)
                mol = Chem.MolFromPDBBlock(pdb_block)
                print('mol 4-1')
                #chem_desc = pypdb.describe_chemical(f"{res_name}")
                print('mol 4-2')
                new_smiles = ''#chem_desc["describeHet"]["ligandInfo"]["ligand"]["smiles"]
                print('mol5')
                template = Chem.MolFromSmiles(new_smiles)
                new_mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
                print('mol6')
                if handle_cov:
                    cov_mol = self.handle_covalent_mol(lig_res_name=res_name, non_cov_mol=new_mol)
                    if cov_mol is not None:
                        new_mol = cov_mol

                return new_mol

            except Exception as e:
                new_pdb_block = ''
                print('mol7')
                for lig in pdb_block.split('\n'):
                    if 'ATM' in lig:
                        pos = 16
                        s = lig[:pos] + ' ' + lig[pos + 1:]
                        new_pdb_block += s
                    else:
                        new_pdb_block += lig

                    new_pdb_block += '\n'
                print('mol8')
                mol = Chem.rdmolfiles.MolFromPDBBlock(new_pdb_block)
                print('mol9')
                return mol

        # Look for new XChem data - new XChem data must have associated smile.txt file
        # Need to do this to catch corner case - x0685 from mArh residue
        # name NHE was found ----> yielded wrong mol/smiles
        if smiles_file:
            new_pdb_block = ''

            for lig in pdb_block.split('\n'):
                if 'ATM' in lig:
                    pos = 16
                    s = lig[:pos] + ' ' + lig[pos + 1:]
                    new_pdb_block += s
                else:
                    new_pdb_block += lig

                new_pdb_block += '\n'

            mol = Chem.rdmolfiles.MolFromPDBBlock(new_pdb_block)
            if handle_cov:
                cov_mol = self.handle_covalent_mol(lig_res_name=res_name, non_cov_mol=mol)
                if cov_mol is not None:
                    mol = cov_mol

            return mol

    def create_pdb_for_ligand(self, ligand, count, monomerize, smiles_file, covalent=False):
        """
        A pdb file is produced for an individual ligand, containing atomic and connection information
        :param ligand: Name of the Ligand
        :param count: The index of the ligand
        :param monomerize: Bool, if the file needs to be named using the chain name of the PDB
        :param smiles_file: File path of smiles_file (if any)
        :param covalent: Bool, indicate whether or not covalent attach should be sought.
        :return: .pdb file for ligand.
        """
        print('lig 1')
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
        print('lig2')
        individual_ligand = []
        individual_ligand_conect = []
        # adding atom information for each specific ligand to a list
        for atom in self.final_hets:
            if str(atom[16:20].strip() + atom[20:26]) == str(ligand):
                individual_ligand.append(atom)

        con_num = 0
        print('lig3')
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
        print('lig4')
        # making into one list that is compatible with conversion to mol object
        ligand_het_con = individual_ligand + individual_ligand_conect
        print('lig5')
        # make a pdb file for the ligand molecule

        if not os.path.isdir(lig_out_dir):
            os.makedirs(lig_out_dir)

        ligands_connections = open(
            os.path.join(lig_out_dir, (file_base + ".pdb")), "w+"
        )
        for line in ligand_het_con:
            ligands_connections.write(str(line))
        ligands_connections.close()
        print('lig6')
        # making pdb file into mol object
        mol = self.create_pdb_mol(file_base=file_base, lig_out_dir=lig_out_dir, smiles_file=smiles_file, handle_cov=covalent)
        print('lig7')
        # Move Map files into lig_out_dir

        if not mol:
            print(f'WARNING: {file_base} did not produce a mol object from its pdb lig file!')
        else:
            try:
                Chem.AddHs(mol)

                self.mol_lst.append(mol)
                self.mol_dict["directory"].append(lig_out_dir)
                self.mol_dict["mol"].append(mol)
                self.mol_dict["file_base"].append(file_base)

            except AssertionError:
                print(file_base, 'is unable to produce a ligand file')
                pass

    def create_mol_file(self, directory, file_base, mol_obj, smiles_file=None):
        """
        a .mol file is produced for an individual ligand
        :param directory: The directory where the mol file should be saved.
        :param file_base: The name of the mol file
        :param mol_obj: The RDKit Mol file object
        :param smiles_file: The filepath of a text file that contains the smiles string of the mol file (if exists).
        :return: A mol file!
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
            print(f'Warning: No smiles file: {file_base}')

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
                smiles = 'NA'

        if not smiles_file:
            try:
                smiles = Chem.MolToSmiles(mol_obj)
            except Exception as e:
                print(e)
                print('failed to convert mol obj to smiles' + smiles_file)
                smiles = "NA"

        meta_data_dict = {'Blank': '',
                          'fragalysis_name': file_base,
                          'crystal_name': file_base.rsplit('_', 1)[0],
                          'smiles': smiles,
                          'new_smiles': '',
                          'alternate_name': '',
                          'site_name': '',
                          'pdb_entry': ''}

        # Write dict to csv
        meta_data_file = open(meta_out_file, 'w+')
        w = csv.DictWriter(meta_data_file, meta_data_dict.keys())
        w.writerow(meta_data_dict)
        meta_data_file.close()


class pdb_apo:
    def __init__(self, infile, target_name, RESULTS_DIRECTORY, filebase, biomol=None):
        self.target_name = target_name
        self.pdbfile = open(infile).readlines()
        self.RESULTS_DIRECTORY = RESULTS_DIRECTORY
        self.filebase = filebase
        self.non_ligs = json.load(
            open(os.path.join(os.path.dirname(__file__), "non_ligs.json"), "r")
        )
        self.apo_file = None
        self.biomol = biomol

    def make_apo_file(self, keep_headers=False):
        """
        Keeps anything other than unique ligands

        :param: pdb file
        :returns: created XXX_apo.pdb file
        """
        lines = ""
        if keep_headers:
            include = ['CONECT', 'SEQRES', 'TITLE', 'ANISOU']
        else:
            include = ['CONECT', 'REMARK', 'CRYST', 'SEQRES', 'HEADER', 'TITLE', 'ANISOU']

        for line in self.pdbfile:
            if (
                    line.startswith("HETATM")
                    and line.split()[3] not in self.non_ligs
                    or any([line.startswith(x) for x in include])
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

        if self.biomol is not None:
            self.add_biomol_remark()
        else:
            print('Not Attaching biomol')

    def add_biomol_remark(self):
        '''
        Add contents of biomol/additional text file to a .pdb file
        '''
        biomol_remark = open(self.biomol).readlines()
        print(biomol_remark)
        f = self.apo_file
        with open(f) as handle:
            switch = 0
            header_front, header_end = [], []
            pdb = []
            for line in handle:
                if line.startswith('ATOM'): switch = 1
                if line.startswith('HETATM'): switch = 2
                if switch == 0:
                    header_front.append(line)
                elif (switch == 2) and not line.startswith('HETATM'):
                    header_end.append(line)
                else:
                    pdb.append(line)
            full_file = ''.join(header_front) + ''.join(biomol_remark) + ''.join(pdb) + ''.join(header_end)
            with open(f, 'w') as w:
                w.write(full_file)

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


def set_up(target_name, infile, out_dir, monomerize, smiles_file=None, biomol=None, covalent=False, keep_headers=False):
    """
    For each ligand inside a pdb file, process each ligand seperately and create own outputs in individual folders.
    :param target_name: Name of the folder in out_dir
    :param infile: pdb file to be processed
    :param out_dir: Overall location of outputs
    :param monomerize: Bool, indicate whether or not data was monomerized prior to set_up
    :param smiles_file: Filepath pointing to text file containing smiles string (if exists)
    :param biomol: Filepath pointing to text file containing biomol/header information for pdbs (if exists)
    :param covalent: Bool, indicate whether or not output mol files should find covalent attachment.
    :param keep_headers: Bool, indicate whether or not keep headers on apo files.
    :return: for each ligand: pdb, mol, sdf and _apo.pdb in seperate directorys inside out_dir/target_name
    """
    print('set1')
    RESULTS_DIRECTORY = os.path.join(out_dir, target_name, 'aligned')

    if not os.path.isdir(RESULTS_DIRECTORY):
        os.makedirs(RESULTS_DIRECTORY)

    print(RESULTS_DIRECTORY)

    new = Ligand(
        target_name, infile, RESULTS_DIRECTORY
    )  # takes in pdb file and returns specific ligand files
    print('set1-1')
    new.hets_and_cons()  # takes only hetatm and conect file lines from pdb file
    print('set1-2')
    new.remove_nonligands()  # removes ions and solvents from list of ligands
    print('set1-3')
    new.find_ligand_names_new()  # finds the specific name and locations of desired ligands
    print('set1-4')
    print(new.wanted_ligs)
    for i in range(len(new.wanted_ligs)):
        new.create_pdb_for_ligand(
            new.wanted_ligs[i], count=i, monomerize=monomerize, smiles_file=smiles_file, covalent=covalent
        )  # creates pdb file and mol object for specific ligand
    print('set2')
    for i in range(len(new.mol_dict["directory"])):
        print('set3')
        if not new.mol_dict["mol"][i]:
            warnings.warn(
                str(
                    "RDkit mol object for "
                    + new.mol_dict["file_base"][i]
                    + " is None, please check the input. Will not write any files"
                )
            )
            continue
        print('set4')
        shutil.copy(infile,
                    os.path.join(new.mol_dict["directory"][i], str(new.mol_dict["file_base"][i] + "_bound.pdb")))

        inpath = infile.replace('_bound.pdb', '')
        basebase = os.path.basename(inpath)
        fofcmap_files = glob.glob(f'{inpath}_*.map')
        event_files = glob.glob(f'{inpath}_*.ccp4')
        map_files = fofcmap_files + event_files
        for map_file in map_files:
            map_base = os.path.basename(map_file)
            map_base = map_base.replace(basebase, new.mol_dict["file_base"][i])
            shutil.copy(map_file,
                        os.path.join(new.mol_dict["directory"][i], map_base))
        print('set5')
        new.create_mol_file(
            directory=new.mol_dict["directory"][i],
            file_base=new.mol_dict["file_base"][i],
            mol_obj=new.mol_dict["mol"][i],
            smiles_file=smiles_file,
        )  # creates mol file for each ligand
        print('set6')
        writer = Chem.rdmolfiles.SDWriter(
            os.path.join(
                new.mol_dict["directory"][i],
                str(new.mol_dict["file_base"][i] + ".sdf"),
            )
        )
        print('set7')
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
        print('set8')
        new_apo = pdb_apo(
            infile,
            target_name,
            new.mol_dict["directory"][i],
            new.mol_dict["file_base"][i],
            biomol=biomol
        )
        new_apo.make_apo_file(keep_headers=keep_headers)  # creates pdb file that doesn't contain any ligand information
        new_apo.make_apo_desol_files()  # makes apo file without solvent, ions and buffers, and file with just those

    return new


def convert_small_AA_chains(in_file, out_file, max_len=15):
    pdb_file = gemmi.read_structure(in_file)
    structure = pdb_file[0]
    chain_lens = [j.get_polymer().length() for j in structure]
    print(chain_lens)
    for i, j in enumerate(chain_lens):
        if int(j) <= int(max_len):
            for x in structure[i]:
                x.name = 'LIG'
                x.het_flag = 'H'
                x.seqid.num = i + 1
    w_chain = structure.find_chain('W')
    nhoh = 1
    if w_chain is not None:
        newchain = gemmi.Chain('W')
        for res in w_chain:
            if res.name == 'LIG' or res.name == 'HOH': # Should be a ligand in the water chain right??
                for atom in res:
                    emptyres = gemmi.Residue()
                    emptyres.add_atom(atom)
                    emptyres.name = 'HOH'
                    emptyres.seqid.num = nhoh
                    nhoh += 1
                    newchain.add_residue(emptyres)
            else:
                newchain.add_residues(res)
        structure.remove_chain('W')
        structure.add_chain(newchain)
    print(out_file)
    pdb_file.write_pdb(out_file)

def copy_extra_files(in_file, out_dir):
    bn = os.path.basename(in_file).replace('.pdb', '')
    dirname = os.path.dirname(in_file)
    other_files = set(glob.glob(os.path.join(dirname, f'{bn}*'))) - set([in_file])
    new_loc = [os.path.join(out_dir, os.path.basename(x)) for x in other_files]
    for i, j in zip(other_files, new_loc):
        shutil.copyfile(i, j)
