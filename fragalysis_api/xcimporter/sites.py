from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors
import numpy as np
import glob
import os
from frag.alysis.run_clustering import run_lig_cluster
import json
import argparse
import dataclasses
import csv

# This needs refactoring.


def get_distance(coord_a, coord_b):
    sum_ = (sum([(float(coord_a[i])-float(coord_b[i]))**2 for i in range(3)]))
    return np.sqrt(sum_)


def json_file_generator(folder):
    jsons = glob.glob(os.path.join(folder, 'aligned', '*', '*_sites.json'))
    if len(jsons) == 0:
        jsons = glob.glob(os.path.join(folder, '*', '*_sites.json'))
    if len(jsons) == 0:
        raise RuntimeError(f'Cannot find any .json files in {folder}')
    return jsons


def mol_file_generator(folder):
    mols = glob.glob(os.path.join(folder, 'aligned', '*', '*.mol'))
    if len(mols) == 0:
        mols = glob.glob(os.path.join(folder, '*', '*.mol'))
    if len(mols) == 0:
        raise RuntimeError(f'Cannot find any .mol files in {folder}')
    return mols


def read_jfile(j_file):
    with open(j_file) as jfile:
        data = json.load(jfile)
    data['fn'] = j_file.replace('_sites.json', '.mol')  # Needed?
    return data


def recluster_from_jsons(folder):
    # Fetch list of all jsons
    jsons = [read_jfile(x) for x in json_file_generator(folder)]
    # Curate dictionary of sites from that
    out_data = {}
    for jf in jsons:
        for keys in jf.keys():
            if keys == 'fn':
                continue
            if out_data.get(keys) is None:
                out_data[keys] = {}
            for n, id in enumerate(jf.get(keys).get('site_id')):
                if out_data[keys].get(id) is None:
                    out_data[keys][id] = {
                        'centre_of_mass': jf.get(keys).get('site_centre_of_mass')[n],
                        'mol_ids': [jf.get('fn')]
                    }
                else:
                    out_data[keys][id]['mol_ids'].append(jf.get('fn'))
    return out_data


def cluster_all_mols(folder):
    mols = mol_file_generator(folder)
    molecules = [Chem.MolFromMolFile(x) for x in mols]
    out_data = run_lig_cluster(molecules, mols)
    return out_data


@dataclasses.dataclass()
class Sites:
    out_data: dict
    mol_files: list
    folder: str
    missing_mols: list

    @staticmethod
    def from_folder(folder, recalculate=True):
        missing_mols = []
        if recalculate:
            outdata = cluster_all_mols(folder)
        else:
            try:
                outdata = recluster_from_jsons(folder)
                missing_mols = [x for x in mol_file_generator(
                    folder) if not os.path.exists(x.replace('.mol', '_sites.json'))]
            except RuntimeError:
                outdata = cluster_all_mols(folder)
        return Sites(out_data=outdata, missing_mols=missing_mols, mol_files=mol_file_generator(folder), folder=folder)

    def to_json(self):
        molids = []
        for key in self.out_data.keys():
            for site in self.out_data.get(key).keys():
                molids = molids + \
                    self.out_data.get(key).get(site).get('mol_ids')
        for x in set(molids):
            out_dict = {}
            for y in self.out_data.keys():
                site_ids = [z for z in self.out_data.get(
                    y).keys() if x in self.out_data.get(y).get(z).get('mol_ids')]
                site_centre_of_mass = [self.out_data.get(y).get(z).get('centre_of_mass') for z in self.out_data.get(
                    y).keys() if x in self.out_data.get(y).get(z).get('mol_ids')]
                if(len(site_ids) > 0):
                    out_dict[y] = {'site_id': site_ids,
                                   'site_centre_of_mass': site_centre_of_mass}
            fn = x.replace('.mol', '_sites.json')
            print(f'Writing Sites to {fn}...')
            with open(fn, 'w') as f:
                json.dump(out_dict, f)

    def apply_to_metadata(self):
        for mf in self.mol_files:
            mf_csv = mf.replace('.mol', '_meta.csv')
            mf_json = mf.replace('.mol', '_sites.json')
            jf = read_jfile(mf_json)
            site_str = str(jf['c_of_m']['site_id'][0] + 1).zfill(2)
            with open(mf_csv, 'r') as csvfile:
                reader = csv.reader(csvfile)
                row = [x for x in reader][0]
            row[6] = f'Site {site_str}'
            with open(mf_csv, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter=',')
                writer.writerow(row)

    def curate_tags(self):
        # Currently undefined, placeholder
        pass

    def cluster_get_or_create(self, key):
        created = False
        if self.out_data.get(key) is None:
            created = True
            self.out_data[key] = {}
        return (self.out_data.get(key), created)

    def next_site_number(self, key):
        vals = list(self.out_data[key].keys())
        site_num = 0
        while site_num in vals:
            site_num += 1
        return site_num

    def tolerate_sites(self, key, com, tolerance, fn):
        cluster, created = self.cluster_get_or_create(key=key)
        if not created:
            if len(cluster.keys()) > 0:
                coms = [cluster.get(x).get('centre_of_mass')
                        for x in cluster.keys()]
                dists = [get_distance(com, y) for y in coms]
                trues = [x < tolerance for x in dists]
                if any(trues):
                    idx = np.argmin(dists)
                    key_idx = list(cluster.keys())[idx]
                    cluster.get(key_idx).get('mol_ids').append(fn)
                else:
                    site_num = self.next_site_number(key=key)
                    cluster[site_num] = {
                        'centre_of_mass': com, 'mol_ids': [fn]}
            else:
                cluster[0] = {'centre_of_mass': com, 'mol_ids': [fn]}
        else:
            cluster[0] = {'centre_of_mass': com, 'mol_ids': [fn]}

    def append_data(self, dict, filename, com_tolerance=5, other_tolerance=1):
        for key in dict.keys():
            if key == 'c_of_m':
                tolerance = com_tolerance
            else:
                tolerance = other_tolerance
            for s in dict.get(key).keys():
                for com in [dict.get(key).get(s).get('centre_of_mass')]:
                    self.tolerate_sites(
                        key=key, com=com, tolerance=tolerance, fn=filename)

    def append_mol(self, mol, com_tolerance=5, other_tolerance=1):
        x_out_data = run_lig_cluster([Chem.MolFromMolFile(mol)], [mol])
        self.append_data(dict=x_out_data, com_tolerance=com_tolerance,
                         other_tolerance=other_tolerance, filename=mol)

    def append_json(self, jsonfile, com_tolerance=5, other_tolerance=1):
        jf = read_jfile(jsonfile)
        for key in jf.keys():
            if key == 'fn':
                continue
            if key == 'c_of_m':
                tolerance = com_tolerance
            else:
                tolerance = other_tolerance
                for com in dict.get(key).get('site_centre_of_mass'):
                    self.tolerate_sites(
                        key=key, com=com, tolerance=tolerance, fn=jf.get('fn'))

    def cluster_missing_mols(self, com_tolerance=5, other_tolerance=1):
        for x in self.missing_mols:
            self.append_mol(mol=x, com_tolerance=com_tolerance,
                            other_tolerance=other_tolerance)


def centre_of_mass(mol):
    numatoms = mol.GetNumAtoms()
    conf = mol.GetConformer()
    if not conf.Is3D():
        return 0
    # get coordinate of each atoms
    pts = np.array([list(conf.GetAtomPosition(atmidx))
                   for atmidx in range(numatoms)])
    atoms = [atom for atom in mol.GetAtoms()]
    mass = Descriptors.MolWt(mol)
    # get center of mass
    center_of_mass = np.array(
        sum(atoms[i].GetMass() * pts[i] for i in range(numatoms))) / mass
    return center_of_mass


def compare_2_mols(molfile, pair):
    # Imagine being this obtuse.
    # , isomericSmiles=False, sanitize=True, canonicalize=True)
    mol1, mol2 = (Chem.MolFromMolFile(x) for x in [molfile, pair])
    smi1, smi2 = (Chem.MolToSmiles(m) for m in [mol1, mol2])
    fps = [Chem.RDKFingerprint(Chem.MolFromSmiles(x)) for x in [smi1, smi2]]
    similarity = DataStructs.FingerprintSimilarity(fps[0], fps[1])
    distance = get_distance(centre_of_mass(mol1), centre_of_mass(mol2))
    num1, num2 = (x.GetNumAtoms() for x in [mol1, mol2])
    conf1, conf2 = (x.GetConformer() for x in [mol1, mol2])
    pts1 = [list(conf1.GetAtomPosition(atmidx)) for atmidx in range(num1)]
    pts2 = [list(conf2.GetAtomPosition(atmidx)) for atmidx in range(num2)]
    dist_list = []
    for p in pts1:
        for q in pts2:
            dist_list.append(get_distance(p, q))
            min_d = min(dist_list)
    # Some rules:
    # If similarity >= .95 + closeatom_Dictance = Really Small) : Alternative Conformation # If same chain Conf, else Pose
    # If similarity >= .95 + closeatom_Distance = Big : Alternative Position
    # If similarity < % + closeatomDistance = big, Primary
    # If similarity < % + closeatomDistance = big, Secondary
    if similarity >= 0.95:
        if min_d <= 5:  # 5 Angstroms enough?
            if molfile.replace('.mol', '')[-1] == pair.replace('.mol', '')[-1]:
                relationship = 'Alternative Conformation'
            else:
                relationship = 'Alternative Pose'
        else:
            relationship = 'Alternative Position'
    else:
        if min_d <= 5:
            relationship = 'Primary'
        else:
            relationship = 'Secondary'
    return (similarity, distance, min_d, relationship)


def contextualize_a_crystal(crystal_list):
    for molfile in crystal_list:
        d = {}
        for pair in list(set(crystal_list) - set([molfile])):
            similarity, distance, min_d, relationship = compare_2_mols(
                molfile, pair)
            d[os.path.basename(pair).replace('.mol', '')] = {
                'similarity': similarity,
                'com_Distance': distance,
                'closeatom_Distance': min_d,
                'relationship': relationship
            }
        fn = molfile.replace('.mol', '_relationships.json')
        with open(fn, 'w') as f:
            json.dump(d, f)


def contextualize_crystal_ligands(folder):
    mols = mol_file_generator(folder)
    unique_mols = set([os.path.basename(x).rsplit('_', 1)[0] for x in mols])
    dic = {}
    for i in unique_mols:
        dic[i] = [x for x in mols if i in x]
    for crys in dic.keys():
        if len(dic.get(crys)) > 1:
            contextualize_a_crystal(dic.get(crys))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--input",
        help="Input folder",
        required=True
    )

    parser.add_argument(
        "-m",
        "--mode",
        choices=['new', 'missing'],
        required=True
    )

    parser.add_argument(
        "-om",
        "--overwritemeta",
        action="store_true",
        help='Include this flag if you would like to automatically assign center of mass sites as site labels',
        required=False,
        default=False
    )

    parser.add_argument(
        "-cc",
        "--cocrystal",
        action="store_true",
        help='Include this flag if you want to define relationships between all ligands within each crystal',
        required=False,
        default=False
    )

    parser.add_argument(
        "-ct",
        "--comtolerance",
        help="Tolerance value for creating new clusters for centre of mass sites",
        type=float,
        default=5.00
    )

    parser.add_argument(
        "-ot",
        "--othertolerance",
        help="Tolerance value for creating new clusters for non centre of mass sites",
        type=float,
        default=1.00
    )

    args = vars(parser.parse_args())
    folder = args["input"]
    mode = args["mode"]
    ot = args["othertolerance"]
    ct = args["comtolerance"]
    cc = args["cocrystal"]
    om = args["overwritemeta"]

    # If recalculate is True then stuff happens, if False it reads from jsons.
    site_obj = Sites.from_folder(folder, recalculate=mode == 'new')
    if mode == 'missing':
        site_obj.cluster_missing_mols(com_tolerance=ct, other_tolerance=ot)
    site_obj.to_json()
    if cc:
        contextualize_crystal_ligands(folder=folder)  # Add arguments?
    if om:
        site_obj.apply_to_metadata()
