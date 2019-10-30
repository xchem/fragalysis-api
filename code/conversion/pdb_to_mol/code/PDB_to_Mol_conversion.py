#!/usr/bin/env python

from rdkit import Chem
import json

pdbcode = "5qj4"
pdbfile = open("../data/"+str(pdbcode)+".pdb").readlines()
non_ligs = json.load(open("non_ligs", "r"))
writer = Chem.rdmolfiles.SDWriter("../resulting_data/"+str(pdbcode)+"_out.sdf")


def hets_and_cons(pdbfile):
    """
    Heteroatoms and connect files are pulled out from full pdb file
    
    param: pdb file in .readlines() format
    returns: lits of hetatomic information and connection information
    """
    hetatms  = []
    conects = []

    for line in pdbfile:
        if line.startswith("HETATM"):
            hetatms.append(line)
        if line.startswith("CONECT"):
            conects.append(line)
            
    return hetatms, conects

           
def remove_nonligands(hetatms, non_ligs):
    """
    Non-ligands such as solvents and ions are removed from the list of heteroatoms, 
    to ideally leave only the target ligand
    
    params: list of heteroatoms and their information, list of non-ligand small molecules 
        that could be in crystal structure
    retunrns: list of heteroatoms that are not contained in the non_ligs list
    """
    final_hets = []
    for line in hetatms:
        if str(line.split()[3]) not in non_ligs:
            final_hets.append(line)
    return final_hets


def find_ligand_names_new(pdbfile, non_ligs):
    """
    Finds list of ligands contained in the structure, including 
    """
    all_ligands = []
    for line in pdbfile:
        if line.startswith("HET "):
            all_ligands.append(line)
    wanted_ligs = []
    for lig in all_ligands:
        if lig.split()[1] not in non_ligs:
            wanted_ligs.append(lig.split()[1:-1])
    return wanted_ligs


def create_pdb_for_ligand(pdbcode, ligand, final_hets, conects):
    """
    A pdb file is produced for an individual ligand, containing atomic and connection information
    
    params: vari pdb conversion, ligand definition, list of ligand heteroatoms and information, connection information
    returns: .pdb file for ligand
    """
    if len(ligand) == 3:
        ligand_name = str(ligand[0]+'_'+str(ligand[1])+'_'+str(ligand[2]))
    elif len(ligand) == 2:
        ligand_name = str(ligand[0]+'_'+str(ligand[1]))
    else:
        ligand_name = str(ligand)
    
    ligands_connections = open("../resulting_data/"+str(pdbcode)+"_"+str(ligand_name)+".pdb", "w+")
    individual_ligand = []
    individual_ligand_conect = []
   
    print(ligand)

    if len(ligand) == 2:
        for atom in final_hets:
            if len(atom.split()[2]) <= 3 and ligand[0] in atom.split()[3] and ligand[1] == atom.split()[4]:
                individual_ligand.append(atom)
            elif len(atom.split()[2]) > 3 and ligand[0] in atom.split()[2] and ligand[1] == atom.split()[3]:
                individual_ligand.append(atom)
    if len(ligand) == 3:
        for atom in final_hets:
            if len(atom.split()[2]) <= 3 and ligand[0] in atom.split()[3] and ligand[1] == atom.split()[4] and ligand[2] == atom.split()[5]:
                individual_ligand.append(atom)
            elif len(atom.split()[2]) > 3 and ligand[0] in atom.split()[2] and ligand[1] == atom.split()[3] and ligand[2] == atom.split()[4]:
                individual_ligand.append(atom)

    print(len(individual_ligand))
            
    for atom in individual_ligand:
        atom_number  = atom.split()[1]
        for conection in conects:
            if atom_number in conection and conection not in individual_ligand_conect:
                individual_ligand_conect.append(conection)

    ligand_het_con = individual_ligand+individual_ligand_conect
    for line in ligand_het_con:
        ligands_connections.write(str(line))
    ligands_connections.close()
    
    return Chem.rdmolfiles.MolFromPDBFile("../resulting_data/"+str(pdbcode)+"_"+str(ligand_name)+".pdb")


def create_mol_file(ligand, mol_obj, pdbcode):
    """
    a .mol file is produced for an individual ligand
    
    params: ligand definition, pdb file, pdb conversion
    returns: .mol file for the ligand
    """
    if len(ligand) == 3:
        ligand_name = str(ligand[0]+'_'+str(ligand[1])+'_'+str(ligand[2]))
    elif len(ligand) == 2:
        ligand_name = str(ligand[0]+'_'+str(ligand[1]))
    else:
        ligand_name = str(ligand)

    return Chem.rdmolfiles.MolToMolFile(mol_obj, "../resulting_data/"+str(pdbcode)+"_"+str(ligand_name)+"_mol.mol")


def create_sd_file(mol_obj, writer):
    """
    a molecular object defined in the pdb file is used to produce a .sdf file
    
    params: pdb file for the molecule, SDWriter from rdkit
    returns: .sdf file with all input molecules from each time the function is called
    """
    return writer.write(mol_obj)            


def main():
    hetatms, conects = hets_and_cons(pdbfile)
    final_hets = remove_nonligands(hetatms, non_ligs)
    ligands = find_ligand_names_new(pdbfile, non_ligs)

    mol_list = []
    for i in range(len(ligands)):
        mol_list.append(create_pdb_for_ligand(pdbcode, ligands[i], final_hets, conects))
    for i in range(len(mol_list)):
        create_mol_file(ligands[i], mol_list[i], pdbcode)
        create_sd_file(mol_list[i], writer)
    writer.close()  # this is important to make sure the file overwrites


main()
