import os
import math
from fragalysis_api.xcwonka.cluster import dp_means
from rdkit.Chem import ChemicalFeatures
from rdkit import Chem
import argparse

# Generate ph4s for a molecule
class RDKitPh4(object):

    factory = None

    def __init__(self):
        # Generate the factory on init
        self.get_factory()

    def get_factory(self):
        """
        Generate the Ph4 feature factory
        :return:
        """
        if self.factory is None:
            this_dir, this_filename = os.path.split(__file__)
            data_path = os.path.join(this_dir, "RDKitPh4.fdef")
            self.factory = ChemicalFeatures.BuildFeatureFactory(data_path)
        return self.factory

    def generate_ph4_for_mol(self, rdmol):
        """
        Generate a pharmacophore from an input molecule and a feature factory.
        :param rdmol: the input RDKit molecule
        :param factory: the feature factory
        :return: a list of 4 tuples (x,y,z, feature)
        """
        feats = self.get_factory().GetFeaturesForMol(rdmol)
        return [
            (feat.GetPos().x, feat.GetPos().y, feat.GetPos().z, feat.GetType())
            for feat in feats
        ]


class RDKitAtom(object):

    def generate_atoms_for_mol(self, rdmol):
        """
        Generate the atoms from an input molecule and a feature factory.
        :param rdmol: the input RDKit molecule
        :return: a list of 4 tuples (x,y,z, atom description)
        """
        out_list = []
        conf = rdmol.GetConformer()
        for atom in rdmol.GetAtoms():
            atom_desc = self.get_atom_description(atom)
            atom_pos = conf.GetAtomPosition(atom.GetIdx())
            out_list.append((atom_pos.x, atom_pos.y, atom_pos.z, atom_desc))
        return out_list

    def get_atom_description(self, atom):
        """
        Generate a unique description of an atom
        :param atom: the input atom
        :return: a hash string of atomic number and  hybridization state
        """
        return "_".join(
            [str(x) for x in [atom.GetAtomicNum(), atom.GetHybridization()]]
        )


def _get_c_of_mass(rdmol):
    """
    Get the unweighted centre of mass of an RDKit Molecule
    :param rdmol:
    :return:
    """
    atoms = rdmol.GetAtoms()
    conf = rdmol.GetConformer()
    x_coord = y_coord = z_coord = 0.0
    numatoms = 0.0
    # Assume all heavy atoms have the same mass
    for atom in atoms:
        if atom.GetAtomicNum() == 1 or atom.GetSmarts() == "[*]":
            continue
        numatoms += 1.0
        coords = conf.GetAtomPosition(atom.GetIdx())

        x_coord += float(coords.x)
        y_coord += float(coords.y)
        z_coord += float(coords.z)
    # Now we have all the coords -> we want to loop through
    if numatoms == 0:
        raise ValueError("No atoms in Molecules")
    return x_coord / numatoms, y_coord / numatoms, z_coord / numatoms


def _parse_pdb(data):
    return Chem.MolFromPDBFile(data)

#PH4_LAMBDA = 1.0
PH4_LAMBDA = 2.0
C_OF_M_LAMBDA = 6.0


def parse_ligand_ph4s(input_mols):
    """
    Function to return a series of ligand based pharmacophores.
    :param input_mols: the RDKit molecules
    :return: the molecule based pharmacophores
    """
    rdkit_ph4 = RDKitPh4()
    rdkit_atom = RDKitAtom()
    output_pharma_list = []
    for mol in input_mols:
        if not mol:
            pharma_list = []
        else:
            pharma_list = rdkit_ph4.generate_ph4_for_mol(rdmol=mol)
            atom_list = rdkit_atom.generate_atoms_for_mol(mol)
            x, y, z = _get_c_of_mass(rdmol=mol)
            c_of_m_feat = (x, y, z, "c_of_m")
            pharma_list.append(c_of_m_feat)
            pharma_list.extend(atom_list)
            print(mol, 'added to output ph4 list')
        output_pharma_list.append(pharma_list)
    return output_pharma_list

def build_type_dict(mol_ph4_list, identifiers):
    type_dict = {}
    for i, mol in enumerate(mol_ph4_list):
        for ph4 in mol:
            x = ph4[0]
            y = ph4[1]
            z = ph4[2]
            ph4_type = ph4[3]
            if ph4_type in type_dict:
                type_dict[ph4_type]["coords"].append((x, y, z))
                type_dict[ph4_type]["mols"].append(identifiers[i])
            else:
                type_dict[ph4_type] = {"coords": [(x, y, z)], "mols": [identifiers[i]]}
    return type_dict


def map_cluster(dp_means_cluster, mol_id_list):
    """

    :param dp_means_cluster:
    :param mol_id_list:
    :return:
    """
    out_dict = {}
    for cluster in dp_means_cluster.clusters:
        out_dict[cluster] = {
            "centre_of_mass": dp_means_cluster.clusters[cluster],
            "mol_ids": [],
        }
    for i, cluster_id in enumerate(dp_means_cluster.dataClusterId):
        out_dict[cluster_id]["mol_ids"].append(mol_id_list[i])
    return out_dict


def cluster_dp(vect_list, lam, mol_list):
    """
    Perform a DP Means clustering.
    :param vect_list: a list of lists of vectors
    :param lam: the clustering parameters
    :param mol_list: the molecular identifers in the same order as the list of vectors
    :return: a dictionary of the form {cluster_id: {centre_of_mass: (x,y,z), mol_ids: [1,5,12]}}
    """
    return map_cluster(dp_means(vect_list, lam), mol_list)


def run_lig_cluster(mols, identifiers):
    """
    Cluster a list of molecules
    :param mols: the input list of RDKit molecules
    :param identifiers: the corresponding list of identifiers
    :return:
    """
    # First we get the list of mols with their Ph4s
    print('parsing ligand ph4s')
    mol_ph4_list = parse_ligand_ph4s(mols)
    print('created ph4s')
    # Then we build a dict of type: coords: [coords list], mols: [mol_index]
    type_dict = build_type_dict(mol_ph4_list, identifiers)
    print('built type dictionary')
    # Then we cluster coords
    clusters = {}
    for ph4_type in type_dict:
        if ph4_type == "c_of_m":
            clusters[ph4_type] = cluster_dp(
                type_dict[ph4_type]["coords"],
                C_OF_M_LAMBDA,
                type_dict[ph4_type]["mols"],
            )
        else:
            clusters[ph4_type] = cluster_dp(
                type_dict[ph4_type]["coords"], PH4_LAMBDA, type_dict[ph4_type]["mols"]
            )
    print('clustered')
    return clusters


def lig_cluster(DATA_DIRECTORY: str = os.path.abspath('data')):
    '''
    clusters by ligands
    :param DATA_DIRECTORY: the directory in which the data director is in
    :return output: clusters of ligands 
    '''
   
    print('directory: ', DATA_DIRECTORY)

    for dir in os.listdir(DATA_DIRECTORY):

        print('directory: ', dir)
        mols = []
        identifiers = []

        for file in os.listdir(os.path.join(DATA_DIRECTORY, dir)):

            if 'bound' in file:

                print('bound_file: ', file)
                mols.append(_parse_pdb(os.path.join(DATA_DIRECTORY, dir, file)))
                identifiers.append(str(file))

        output = run_lig_cluster(mols, identifiers)

    return output

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-id', '--user_id', required=True,
                        help='Description for foo argument')
    args = vars(parser.parse_args())

    user_id = args['user_id']

    in_dir =  os.path.join('..', '..', 'data', 'xcimporter', 'output', user_id)

    trial = lig_cluster(in_dir)

    print('final output: ', trial)
    print('The ligands have been clustered!!!!')