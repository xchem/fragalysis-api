import json
import os
from rdkit import Chem
import numpy as np
from rdkit.Chem import AllChem

def get_3d_distance(coord_a, coord_b):
    sum_ = (sum([(float(coord_a[i])-float(coord_b[i]))**2 for i in range(3)]))
    return np.sqrt(sum_)


links = json.load(open('file-smiles.json', 'r'))
target = 'Mpro'
lig = 'x0689'
smile_ref = target+'-'+lig+'-bound.pdb'
pdb = open(os.path.join('data', 'Mpro_allPdb_01-Apr-2020', target+'-'+lig+'_0.pdb'), 'r').readlines()
f = open('Mpro_lig_2.pdb', 'w')
link_atoms = {'SG':16, 'O': 8, 'N': 7}

lig_atoms = []

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

m = Chem.MolFromSmiles(links[smile_ref])
m2 = Chem.MolFromPDBFile('Mpro_lig_2.pdb')
print(m, m2)
print(link_info)

for line in pdb:
    if 'CYS' in line and line.startswith('ATOM'):
        if line[13:17].strip() == link_info[1] and line[17:20] == link_info[2] and line[20:23].strip() == link_info[3] and line[23:27].strip() == link_info[4]:
            res_coords = [line[31:39].strip(), line[39:47].strip(), line[47:55].strip()]


edmol = Chem.EditableMol(m2)
new_mol = edmol.AddAtom(Chem.Atom(link_atoms[link_info[1]]))
Chem.MolToPDBFile(edmol.GetMol(), 'edlig.pdb')
edpdb = open('edlig.pdb', 'r').readlines()

distances = {}
for line in edpdb:
    if line.startswith('HETATM'):
        distances[get_3d_distance(res_coords, [line[31:39].strip(), line[39:47].strip(), line[47:55].strip()])] = line[7:11].strip()
    if [line[31:39].strip(), line[39:47].strip(), line[47:55].strip()] == ['0.000', '0.000', '0.000']:
        new_idx = line[7:11].strip()

print(distances[min(distances)], new_idx)
edmol.AddBond(int(distances[min(distances)])-1, int(new_idx)-1, Chem.BondType.SINGLE)

from rdkit.Geometry import Point3D
new_mol = edmol.GetMol()
conf = new_mol.GetConformer()
res_coords = tuple([float(i) for i in res_coords])
conf.SetAtomPosition(int(new_idx)-1, Point3D(res_coords[0], res_coords[1], res_coords[2]))
print(Chem.MolToPDBBlock(new_mol))

