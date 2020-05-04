import json
import os
from rdkit import Chem
import numpy as np
from rdkit.Chem import AllChem

def get_3d_distance(coord_a, coord_b):
    sum_ = (sum([(float(coord_a[i])-float(coord_b[i]))**2 for i in range(3)]))
    return np.sqrt(sum_)

def write_lig_pdb(pdb):
    f = open('lig.pdb', 'w')
    lig_atoms = []
    link_info = None
    for line in pdb:
        if line.startswith('LINK'):
            f.write(line)
            link_info = [line[13]]+[line[42:46].strip()]+line[47:62].strip().split(' ')
        if line[17:20] == 'LIG' and line.startswith('HETATM'):
            lig_atoms.append(line[7:11].strip())
            f.write(line)
        if line.startswith('CONECT'):
            yes = False
            for num in line.split(' '):
                if num in lig_atoms:
                    yes = True
            if yes == True:
                f.write(line)
    f.close()
    return link_info


target = 'Mpro'
if os.path.exists(os.path.join('results', target)):
    os.system('rm -r '+os.path.join('results', target))
os.makedirs(os.path.join('results', target))

link_atoms = {'SG':16, 'O': 8, 'N': 7, 'S' : 16}

for file in os.listdir(os.path.join('data', target)):
    lig = file[5:10]
    smile_ref = target+'-'+lig+'-bound.pdb'
    pdb = open(os.path.join('data', 'Mpro', target+'-'+lig+'_0.pdb'), 'r').readlines()

    link_info = write_lig_pdb(pdb)
    if link_info is not None:
        m2 = Chem.MolFromPDBFile('lig.pdb')

        for line in pdb:
            if line[13:17].strip() == link_info[1] and line[17:20] == link_info[2] and line[20:23].strip() == link_info[3] and line[23:27].strip() == link_info[4]:
                try:
                    res_coords = [float(line[31:39].strip()), float(line[39:47].strip()), float(line[47:55].strip())]
                except ValueError:
                    pass

        edmol = Chem.EditableMol(m2)
        try:
            new_mol = edmol.AddAtom(Chem.Atom(link_atoms[link_info[1]]))
        except ValueError:
            new_mol = edmol.AddAtom(Chem.Atom(link_atoms[link_info[0]]))
        Chem.MolToPDBFile(edmol.GetMol(), 'edlig.pdb')
        edpdb = open('edlig.pdb', 'r').readlines()

        distances = {}
        for line in edpdb:
            if line.startswith('HETATM'):
               # print(res_coords)
                distances[get_3d_distance(res_coords, [line[31:39].strip(), line[39:47].strip(), line[47:55].strip()])] = line[7:11].strip()
            if [line[31:39].strip(), line[39:47].strip(), line[47:55].strip()] == ['0.000', '0.000', '0.000']:
                new_idx = line[7:11].strip()

        edmol.AddBond(int(distances[min(distances)])-1, int(new_idx)-1, Chem.BondType.SINGLE)

        from rdkit.Geometry import Point3D
        new_mol = edmol.GetMol()
        conf = new_mol.GetConformer()
        res_coords = tuple([float(i) for i in res_coords])
        conf.SetAtomPosition(int(new_idx)-1, Point3D(res_coords[0], res_coords[1], res_coords[2]))
        Chem.MolToMolFile(new_mol, os.path.join('results', target, file[:-3]+'mol'))
