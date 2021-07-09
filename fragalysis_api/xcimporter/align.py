import glob
import time

import Bio.PDB as bp
from pathlib import Path
from scipy import spatial
from pandda_gemmi.pandda_types import *
import dataclasses
import pandas as pd
import os
import warnings
import Bio.PDB.PDBExceptions as bpp
import shutil
import gemmi  # Oh boy...
import numpy as np

warnings.simplefilter('ignore', bpp.PDBConstructionWarning)


class Align:

    def __init__(self, directory, pdb_ref='', rrf=False, refset=True):
        '''
        :param directory: Directory path contain pdbs to be aligned.
        :param pdb_ref: The String Reference of a pdb that you want to
            optionally align the input pdbs to. This should be the name of the PDB without
            .pbd and should be located inside the specified directory.
        :param rrf: Bool, To indicate if chains have been split into seperate pdb files.
            If True, this will change the behaviour of the alignment to solely align the first
            chain of each structure together, and not consider other chains.
        '''

        self.directory = directory
        if refset:
            self._get_ref = pdb_ref
        self.rrf = rrf
        self.ref_map = ''

    @property
    def _get_files(self):
        """
        Extracts a list of filepaths for all PDB files within the input directory
        :return: list of .pdb file names in directory
        """
        return glob.glob(os.path.join(self.directory, '*.pdb'))

    @property
    def _get_maplist(self):
        """
        Extracts a list of all files that get with .ccp4 or .map
        :return: a list of .map or .ccp4 file names in a directory.
        """
        map_files = glob.glob(os.path.join(self.directory, '*.map'))
        ccp4_files = glob.glob(os.path.join(self.directory, '*.ccp4'))
        return map_files + ccp4_files

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
        :return: Sets the name of the __pdb_ref attribute to the align class depending on which pdb was chosen
        """
        if pdb_ref != '':
            try:
                assert (os.path.isfile(os.path.join(
                    self.directory, pdb_ref + ".pdb")) == 1)
                self.__pdb_ref = pdb_ref
            except AssertionError:
                print('pdb desired as reference does not exists. Default pdb chosen.')
                self.__pdb_ref = self.__best_length_and_resolution(
                    [i for i in self._get_files if 'pdb' in i])
        else:
            self.__pdb_ref = self.__best_length_and_resolution(
                [i for i in self._get_files if 'pdb' in i])

    def __best_length_and_resolution(self, pdb_files):
        """
        Find the longest pdb structure with the lowest resolution from all imported files.
        This will be used as the reference pdb for alignment against.
        :param pdb_files: List of pdb filepaths
        :return: str of filename with best .pdb file
        """
        a_df = pd.DataFrame(pdb_files, columns=['file'])
        a_df[['resolution', 'p_len', 'pdb']] = a_df.file.apply(
            self.__get_length_and_resolution)

        return a_df.sort_values(by=['p_len', 'resolution'], ascending=[False, True]).pdb[0]

    def __get_length_and_resolution(self, file):
        """
        Determine resolution, sequence and length of .pdb file.
        :param file: pdb file path.
        :return: pandas series with resolution, sequence and length and .pdb filename.
        """
        parser = bp.PDBParser()
        ppb = bp.PPBuilder()
        structure = parser.get_structure(
            os.path.splitext(os.path.basename(file))[0], file)

        seq_len = 0

        # Retrieve length by looping through each chain in the protein
        for pp in ppb.build_peptides(structure):
            seq_len += len(pp.get_sequence())

        # using a functions from PDBParser parser class to get the resolution and protein id from the pdb file
        return pd.Series([structure.header['resolution'], seq_len, structure.id])

    def align_to_reference(self, in_file, reference_pdb, out_dir):
        '''
        Aligns a single pdb file to a reference and adds it to the specified output directory.
        :param in_file: filepath to corresponding pdb file to align. Accompanying map files should be located within
            the same directory as input file.
        :param reference_pdb: the reference pdb file to align to.
        :param out_dir: The desired output directory for the aligned pdb file.
        :return: an aligned pdb file in the output directory with the same name as input.
        '''
        input_files = in_file
        map_list = self._get_maplist
        base_names = os.path.splitext(os.path.basename(input_files))[0]
        crystals = [y for y in [x for x in [base_names]
                                if 'event' not in x] if 'fofc' not in y]
        ref = reference_pdb
        dir = self.directory
        rrf = self.rrf

        reference_pdb = Structure.from_file(file=Path(ref))
        s = time.time()
        for num, name in enumerate(crystals):
            all_maps = [j for j in map_list if name in j]
            current_pdb = Structure.from_file(file=Path(in_file))
            if rrf:
                # current_pdb.structure, chains = split_chain_str(
                chains = split_chain_str(os.path.join(dir, f'{name}.pdb'))
            else:
                chains = ['']
            for chain in chains:
                if rrf:
                    print(
                        f'Aligning Chain {chain} of {name} to first chain of {ref}')
                try:
                    current_pdb, transform = current_pdb.align_to(
                        other=reference_pdb, rrf=rrf, chain_id=chain
                    )
                except Exception as e:
                    print(f'{e}')
                    continue
                # Write new structure according to chain-name?
                if rrf:
                    current_pdb.structure.write_pdb(
                        os.path.join(out_dir, f'{name}_{chain}_bound.pdb')
                    )
                    if os.path.exists(os.path.join(self.directory, f'{name}_smiles.txt')):
                        shutil.copyfile(os.path.join(self.directory, f'{name}_smiles.txt'), os.path.join(
                            out_dir, f'{name}_{chain}_smiles.txt'))
                else:
                    current_pdb.structure.write_pdb(
                        os.path.join(out_dir, f'{name}_bound.pdb')
                    )
                    if os.path.exists(os.path.join(self.directory, f'{name}_smiles.txt')):
                        shutil.copyfile(os.path.join(self.directory, f'{name}_smiles.txt'), os.path.join(
                            out_dir, f'{name}_smiles.txt'))

                # Align Xmaps + save!
                for i in all_maps:
                    base, ext = os.path.splitext(os.path.basename(i))
                    s2 = time.time()
                    map = Xmap.from_file(
                        file=Path(os.path.join(dir, f'{base}{ext}')))
                    array = np.array(map.xmap, copy=False)
                    array[~np.isfinite(array)] = 0
                    print(transform)
                    print(transform.transform.vec.tolist())
                    print(transform.transform.mat.tolist())
                    newmap = resample(
                        moving_xmap=map, transform=transform, reference_structure=reference_pdb)
                    template = Path(
                        os.path.join(dir, f'{base}{ext}'))
                    if rrf:
                        base = base.replace(name, f'{name}_{chain}')
                    fn = f'{base}{ext}'
                    referenceSave(
                        template_map_path=template,
                        xmap=newmap,
                        path_to_save=Path(os.path.join(out_dir, fn))
                    )
                    e2 = time.time()
                    print(f'{int(e2 - s2)} seconds to transform map...')
                    e = time.time()
                    print(f'Total Running time: {int(e - s) / 60} minutes.')

    def write_align_ref(self, output):
        '''
        Copy the reference pdb structure used for alignment to a new location
        :param output: The corresponding filename to write the reference pdb to
        :return: the pdbfile that was chosen as reference copied to the located specified by output
        '''
        fn = os.path.join(self.directory, f'{self._get_ref}.pdb')
        shutil.copyfile(fn, os.path.join(output, 'reference.pdb'))

    def align(self, out_dir):
        """
        Aligns all pdb structures and map files (if any) using gemmi, and save them to a new directory.
        :param out_dir: directory to save aligned pdbs in
        :return: saves the pdbs + transforms map files (if any!)
        """
        # load ref
        input_files = self._get_files
        map_list = self._get_maplist
        crystals = [os.path.splitext(os.path.basename(f))[0]
                    for f in input_files]
        ref = self._get_ref
        dir = self.directory
        rrf = self.rrf

        # Reference stuff
        reference_pdb = Structure.from_file(
            file=Path(os.path.join(dir, f'{ref}.pdb')))

        s = time.time()
        for num, name in enumerate(crystals):
            all_maps = [j for j in map_list if name in j]
            # Logic to do...
            # Align Chain N to First Chain in Reference
            current_pdb = Structure.from_file(
                file=Path(os.path.join(dir, f'{name}.pdb'))
            )
            if rrf:
                chains = split_chain_str(
                    os.path.join(dir, f'{name}.pdb'))
            else:
                chains = ['']
            for chain in chains:
                try:
                    current_pdb, transform = current_pdb.align_to(
                        other=reference_pdb, rrf=rrf, chain_id=chain
                    )
                except Exception as e:
                    print(f'{e}')
                    continue

                # Write new structure according to chain-name?
                if rrf:
                    current_pdb.structure.write_pdb(
                        os.path.join(out_dir, f'{name}_{chain}_bound.pdb')
                    )
                    if os.path.exists(os.path.join(self.directory, f'{name}_smiles.txt')):
                        shutil.copyfile(os.path.join(self.directory, f'{name}_smiles.txt'), os.path.join(
                            out_dir, f'{name}_{chain}_smiles.txt'))
                else:
                    current_pdb.structure.write_pdb(
                        os.path.join(out_dir, f'{name}_bound.pdb')
                    )
                    if os.path.exists(os.path.join(self.directory, f'{name}_smiles.txt')):
                        shutil.copyfile(os.path.join(self.directory, f'{name}_smiles.txt'), os.path.join(
                            out_dir, f'{name}_smiles.txt'))

                # Align Xmaps + save!
                for i in all_maps:
                    base, ext = os.path.splitext(os.path.basename(i))
                    s2 = time.time()
                    map = Xmap.from_file(
                        file=Path(os.path.join(dir, f'{base}{ext}')))
                    array = np.array(map.xmap, copy=False)
                    array[~np.isfinite(array)] = 0
                    newmap = resample(
                        moving_xmap=map, transform=transform, reference_structure=reference_pdb)
                    template = Path(
                        os.path.join(dir, f'{base}{ext}'))
                    if rrf:
                        base = base.replace(name, f'{name}_{chain}')
                    fn = f'{base}{ext}'
                    print(fn)
                    # newmap.save(
                    #    path=Path(os.path.join(out_dir, fn)))
                    referenceSave(
                        template_map_path=template,
                        xmap=newmap,
                        path_to_save=Path(os.path.join(out_dir, fn))
                    )
                    e2 = time.time()
                    print(f'{int(e2 - s2)} seconds to transform map...')
                    e = time.time()
                    print(f'Total Running time: {int(e - s) / 60} minutes.')


@dataclasses.dataclass()
class ResidueID:
    model: str
    chain: str
    insertion: str

    @staticmethod
    def from_residue_chain(model: gemmi.Model, chain: gemmi.Chain, res: gemmi.Residue):
        return ResidueID(model.name, chain.name, str(res.seqid.num))

    def __hash__(self):
        return hash((self.model, self.chain, self.insertion))


@dataclasses.dataclass()
class RFree:
    rfree: float

    @staticmethod
    def from_structure(structure):
        rfree = structure.structure.make_mmcif_document(
        )[0].find_loop("_refine.ls_R_factor_R_free")[0]
        return RFree(float(rfree))

    def to_float(self):
        return self.rfree


@dataclasses.dataclass()
class Transform:
    transform: gemmi.Transform
    com_reference: np.array
    com_moving: np.array

    def apply(self, position: gemmi.Position) -> gemmi.Position:
        rotation_frame_position = gemmi.Position(position[0] - self.com_reference[0],
                                                 position[1] -
                                                 self.com_reference[1],
                                                 position[2] - self.com_reference[2])
        transformed_vector = self.transform.apply(rotation_frame_position)
        transformed_position = gemmi.Position(transformed_vector[0] + self.com_moving[0],
                                              transformed_vector[1] +
                                              self.com_reference[1],
                                              transformed_vector[2] + self.com_reference[2])
        return transformed_position

    def apply_inverse(self, position: gemmi.Position) -> gemmi.Position:
        rotation_frame_position = gemmi.Position(position[0] - self.com_moving[0],
                                                 position[1] -
                                                 self.com_moving[1],
                                                 position[2] - self.com_moving[2])
        transformed_vector = self.transform.inverse().apply(rotation_frame_position)
        transformed_position = gemmi.Position(transformed_vector[0] + self.com_reference[0],
                                              transformed_vector[1] +
                                              self.com_reference[1],
                                              transformed_vector[2] + self.com_reference[2])
        return transformed_position

    @staticmethod
    def from_translation_rotation(translation, rotation, com_reference, com_moving):
        transform = gemmi.Transform()
        transform.vec.fromlist(translation.tolist())
        transform.mat.fromlist(rotation.as_matrix().tolist())

        return Transform(transform, com_reference, com_moving)

    @staticmethod
    def pos_to_list(pos: gemmi.Position):
        return [pos[0], pos[1], pos[2]]


@dataclasses.dataclass()
class Structure:
    structure: gemmi.Structure

    @staticmethod
    def from_file(file):
        structure = gemmi.read_structure(str(file))
        return Structure(structure)

    def rfree(self):
        return RFree.from_structure(self)

    def __getitem__(self, item: ResidueID):
        return self.structure[item.model][item.chain][item.insertion]

    def residue_ids(self):
        residue_ids = []
        for model in self.structure:
            for chain in model:
                for residue in chain.get_polymer():
                    resid = ResidueID.from_residue_chain(model, chain, residue)
                    residue_ids.append(resid)

        return residue_ids

    def protein_atoms(self):
        for model in self.structure:
            for chain in model:
                for residue in chain.get_polymer():
                    for atom in residue:
                        yield atom

    def all_atoms(self):
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        yield atom

    def align_to(self, other, rrf=False, chain_id=''):
        # TODO: CHANGE
        # Warning: inplace!
        # Aligns structures usings carbon alphas and transform self into the frame of the other

        ca_self = []
        ca_other = []
        ca_res_self_names = []
        ca_res_other_names = []
        # Get CAs
        for model in self.structure:
            for chain in model:
                if rrf and chain.name not in chain_id:
                    continue
                else:
                    for res_self in chain.get_polymer():
                        if 'LIG' in str(res_self):
                            continue
                        try:
                            current_res_id = ResidueID.from_residue_chain(
                                model, chain, res_self)
                            if rrf:  # TODO CHANGE?
                                res_other = other.structure[current_res_id.model][0][current_res_id.insertion][0]
                            else:
                                res_other = \
                                    other.structure[current_res_id.model][current_res_id.chain][current_res_id.insertion][0]
                            self_ca_pos = res_self["CA"][0].pos
                            other_ca_pos = res_other["CA"][0].pos
                            ca_res_self_names.append(res_self)
                            ca_res_other_names.append(res_other)
                        except:
                            continue

                        ca_list_self = Transform.pos_to_list(self_ca_pos)
                        ca_list_other = Transform.pos_to_list(other_ca_pos)

                        ca_self.append(ca_list_self)
                        ca_other.append(ca_list_other)

        print(f'ca_self {len(ca_self)}')
        print(ca_res_self_names)
        print(f'ca_other {len(ca_other)}')
        print(ca_res_other_names)
        # Make coord matricies
        matrix_self = np.array(ca_self)
        matrix_other = np.array(ca_other)

        # Find means
        mean_self = np.mean(matrix_self, axis=0)
        mean_other = np.mean(matrix_other, axis=0)

        # demaen
        de_meaned_self = matrix_self - mean_self
        de_meaned_other = matrix_other - mean_other

        # Align
        rotation, rmsd = spatial.transform.Rotation.align_vectors(
            de_meaned_self, de_meaned_other)

        # Get transform
        vec = np.array([0.0, 0.0, 0.0])
        # Transform is from other frame to self frame
        transform = Transform.from_translation_rotation(vec,
                                                        rotation,
                                                        mean_other,
                                                        mean_self,
                                                        )
        # Transform positions
        for atom in self.all_atoms():
            atom.pos = transform.apply_inverse(atom.pos)

        return self, transform


def split_chain_str(f):
    aa_codes = {'V': 'VAL', 'I': 'ILE', 'L': 'LEU', 'E': 'GLU', 'Q': 'GLN', 'D': 'ASP', 'N': 'ASN', 'H': 'HIS',
                'W': 'TRP', 'F': 'PHE', 'Y': 'TYR', 'R': 'ARG', 'K': 'LYS', 'S': 'SER', 'T': 'THR', 'M': 'MET',
                'A': 'ALA', 'G': 'GLY', 'P': 'PRO', 'C': 'CYS'}
    base_structure = gemmi.read_structure(f)
    base_models = base_structure[0]
    all_chains = [x.name for x in base_models]
    HOH_chains = find_water_chains(f)
    nonHOH_chains = list(set(all_chains) - set(HOH_chains))
    chain_centers = {}
    chain_names = []
    for i in nonHOH_chains:
        chain_centers[i] = get_chain_center(chain_name=i, file=f)
        span = base_structure[0][i].whole()
        residues = [x.name for x in span]
        # Possible to convert this to a %age
        if any([True for v in residues if v in aa_codes.values()]):
            chain_names.append(i)
    alt_chains = list(set(nonHOH_chains) - set(chain_names))
    return list(set(nonHOH_chains) - set(alt_chains))


def split_chains(f):
    '''
    Split a pdb file according to chain names. While preserving chains that contain waters and assigned small, non-amino acid containing chains to the nearest chain by center of mass.
    :param f: File name of the pdb file to split
    :return: A list of pdb_files that have been created with the naming convention f_[chain_name]
    '''
    aa_codes = {'V': 'VAL', 'I': 'ILE', 'L': 'LEU', 'E': 'GLU', 'Q': 'GLN', 'D': 'ASP', 'N': 'ASN', 'H': 'HIS',
                'W': 'TRP', 'F': 'PHE', 'Y': 'TYR', 'R': 'ARG', 'K': 'LYS', 'S': 'SER', 'T': 'THR', 'M': 'MET',
                'A': 'ALA', 'G': 'GLY', 'P': 'PRO', 'C': 'CYS'}
    filenames = []
    base_structure = gemmi.read_structure(f)
    base_models = base_structure[0]
    all_chains = [x.name for x in base_models]
    HOH_chains = find_water_chains(f)
    nonHOH_chains = list(set(all_chains) - set(HOH_chains))

    chain_centers = {}
    chain_names = []
    # From nonHOH chains, identify if they contain any aminoacids, if they do. Assume they are full chains.
    # else, later on they will be attached to nearest chain by center.
    for i in nonHOH_chains:
        chain_centers[i] = get_chain_center(chain_name=i, file=f)
        span = base_structure[0][i].whole()
        residues = [x.name for x in span]
        # Possible to convert this to a %age
        if any([True for v in residues if v in aa_codes.values()]):
            chain_names.append(i)

    alt_chains = list(set(nonHOH_chains) - set(chain_names))
    for i in chain_names:
        # For each chain, convert all ligands,
        temp_structure = gemmi.read_structure(f)
        for j in alt_chains:
            chain_dists = {}
            for z in chain_names:
                chain_dists[z] = chain_centers[j].dist(chain_centers[z])
            temp_structure[0][j].name = min(chain_dists, key=chain_dists.get)
        leftover_chains = list(set(nonHOH_chains) - set(i))
        # Remove remaining chains
        for j in leftover_chains:
            if j == i:
                # Just in case.
                continue
            else:
                temp_structure[0].remove_chain(j)
                #print([x.name for x in temp_structure[0]])

        # Rename Chain to corresponding chain then save!
        name = os.path.splitext(os.path.basename(f))[0] + '_' + str(i)
        filename = os.path.join(self.outdir, f'{name}_mono.pdb')
        #print(f'Writing to: {filename}')
        temp_structure.write_pdb(filename)
        filenames.append(filename)

    return filenames


def get_chain_center(chain_name, file):
    '''
    Calculate the center of mass for a particular chain within a particular pdb file
    :param chain_name: Name of the chain
    :param file: Filepath of pdb file.
    :return: The center of mass of the specified chain.
    '''
    structure = gemmi.read_pdb(file)[0]
    names = [x.name for x in structure]
    for i in names:
        if i == chain_name:
            continue
        else:
            structure.remove_chain(i)
    return structure.calculate_center_of_mass()


def find_water_chains(file):
    '''
    Find which chain(s) contain only water molecules
    :param file: A pdb file name.
    :return: A list of chains that contain water
    '''
    structure = gemmi.read_pdb(file)
    old = [x.name for x in structure[0]]
    structure.remove_waters()
    structure.remove_empty_chains()
    new = [x.name for x in structure[0]]
    return list(set(old) - set(new))


def resample(
        moving_xmap: Xmap,
        transform: Transform,
        reference_structure: Structure
):

    # interpolated_grid = gemmi.FloatGrid(
    #    moving_xmap.xmap.nu,
    #    moving_xmap.xmap.nv,
    #    moving_xmap.xmap.nw, )
    # interpolated_grid.set_unit_cell(moving_xmap.xmap.unit_cell)
    #interpolated_grid.spacegroup = moving_xmap.xmap.spacegroup
    interpolated_grid = gridFromTemplate(moving_xmap)
    # points
    mask = gemmi.FloatGrid(moving_xmap.xmap.nu,
                           moving_xmap.xmap.nv,
                           moving_xmap.xmap.nw, )
    mask.set_unit_cell(moving_xmap.xmap.unit_cell)
    mask.spacegroup = gemmi.find_spacegroup_by_name("P 1")

    for model in reference_structure.structure:
        for chain in model:
            for residue in chain.get_polymer():
                for atom in residue:
                    mask.set_points_around(atom.pos, 5.0, 1.0)

    mask_array = np.array(mask)
    mask_indicies = np.hstack([x.reshape((len(x), 1))
                              for x in np.nonzero(mask)])
    fractional_coords = []
    for model in reference_structure.structure:
        for chain in model:
            for residue in chain.get_polymer():
                for atom in residue:
                    fractional = moving_xmap.xmap.unit_cell.fractionalize(
                        atom.pos)
                    fractional_coords.append(
                        [fractional.x, fractional.y, fractional.z])

    fractional_coords_array = np.array(fractional_coords)
    max_coord = np.max(fractional_coords_array, axis=0)
    min_coord = np.min(fractional_coords_array, axis=0)

    min_index = np.floor(
        min_coord * np.array([interpolated_grid.nu, interpolated_grid.nv, interpolated_grid.nw]))
    max_index = np.floor(
        max_coord * np.array([interpolated_grid.nu, interpolated_grid.nv, interpolated_grid.nw]))

    points = itertools.product(range(int(min_index[0]), int(max_index[0])),
                               range(int(min_index[1]), int(max_index[1])),
                               range(int(min_index[2]), int(max_index[2])),
                               )

    # Unpack the points, poitions and transforms
    point_list: List[Tuple[int, int, int]] = []
    position_list: List[Tuple[float, float, float]] = []
    transform_list: List[gemmi.transform] = []
    com_moving_list: List[np.array] = []
    com_reference_list: List[np.array] = []

    transform_rotate_reference_to_moving = transform.transform
    transform_rotate_reference_to_moving.vec.fromlist([0.0, 0.0, 0.0])

    transform_reference_to_centered = gemmi.Transform()
    transform_reference_to_centered.vec.fromlist(
        (-transform.com_reference).tolist())
    transform_reference_to_centered.mat.fromlist(np.eye(3).tolist())

    transform_centered_to_moving = gemmi.Transform()
    transform_centered_to_moving.vec.fromlist(transform.com_moving.tolist())
    transform_centered_to_moving.mat.fromlist(np.eye(3).tolist())

    # indicies to positions
    for point in points:
        position = interpolated_grid.point_to_position(
            interpolated_grid.get_point(point[0], point[1], point[2]))
        # Tranform to origin frame
        trtc = transform_reference_to_centered.apply(position)
        position_origin_reference = gemmi.Position(trtc[0], trtc[1], trtc[2])

        # Rotate
        trrtm = transform_rotate_reference_to_moving.apply(
            position_origin_reference)
        position_origin_moving = gemmi.Position(trrtm[0], trrtm[1], trrtm[2])

        # Transform to moving frame
        tctm = transform_centered_to_moving.apply(position_origin_moving)
        position_moving = gemmi.Position(tctm[0], tctm[1], tctm[2])

        # Interpolate moving map
        interpolated_map_value = moving_xmap.xmap.interpolate_value(
            position_moving)

        # Set original point
        interpolated_grid.set_value(
            point[0], point[1], point[2], interpolated_map_value)

    # interpolated_grid.symmetrize_max()
    interpolated_array = np.array(interpolated_grid)

    interpolated_grid_neg = gridFromTemplate(moving_xmap)
    interpolated_array_neg = np.array(interpolated_grid_neg, copy=False)
    interpolated_array_neg[:, :, :] = -interpolated_array[:, :, :]
    interpolated_grid_neg.symmetrize_max()

    interpolated_grid_pos = gridFromTemplate(moving_xmap)
    interpolated_array_pos = np.array(interpolated_grid_pos, copy=False)
    interpolated_array_pos[:, :, :] = interpolated_array[:, :, :]
    interpolated_grid_pos.symmetrize_max()

    interpolated_grid_sym = gridFromTemplate(moving_xmap)
    interpolated_array_sym = np.array(interpolated_grid_sym, copy=False)
    interpolated_array_sym[:, :, :] = interpolated_array_pos[:,
                                                             :, :] - interpolated_array_neg[:, :, :]

    return Xmap(interpolated_grid_sym)


def resample2(
        moving_xmap: Xmap,
        transform: Transform,
        reference_structure: Structure
):

    # interpolated_grid = gemmi.FloatGrid(
    #    moving_xmap.xmap.nu,
    #    moving_xmap.xmap.nv,
    #    moving_xmap.xmap.nw, )
    # interpolated_grid.set_unit_cell(moving_xmap.xmap.unit_cell)
    #interpolated_grid.spacegroup = moving_xmap.xmap.spacegroup
    interpolated_grid = gridFromTemplate(moving_xmap)
    # points
    mask = gemmi.FloatGrid(moving_xmap.xmap.nu,
                           moving_xmap.xmap.nv,
                           moving_xmap.xmap.nw, )
    mask.set_unit_cell(moving_xmap.xmap.unit_cell)
    mask.spacegroup = gemmi.find_spacegroup_by_name("P 1")

    for model in reference_structure.structure:
        for chain in model:
            for residue in chain.get_polymer():
                for atom in residue:
                    mask.set_points_around(atom.pos, 5.0, 1.0)

    mask_array = np.array(mask)
    mask_indicies = np.hstack([x.reshape((len(x), 1))
                              for x in np.nonzero(mask)])
    fractional_coords = []
    for model in reference_structure.structure:
        for chain in model:
            for residue in chain.get_polymer():
                for atom in residue:
                    fractional = moving_xmap.xmap.unit_cell.fractionalize(
                        atom.pos)
                    fractional_coords.append(
                        [fractional.x, fractional.y, fractional.z])

    fractional_coords_array = np.array(fractional_coords)
    max_coord = np.max(fractional_coords_array, axis=0)
    min_coord = np.min(fractional_coords_array, axis=0)

    min_index = np.floor(
        min_coord * np.array([interpolated_grid.nu, interpolated_grid.nv, interpolated_grid.nw]))
    max_index = np.floor(
        max_coord * np.array([interpolated_grid.nu, interpolated_grid.nv, interpolated_grid.nw]))

    points = itertools.product(range(int(min_index[0]), int(max_index[0])),
                               range(int(min_index[1]), int(max_index[1])),
                               range(int(min_index[2]), int(max_index[2])),
                               )

    # Unpack the points, poitions and transforms
    point_list: List[Tuple[int, int, int]] = []
    position_list: List[Tuple[float, float, float]] = []
    transform_list: List[gemmi.transform] = []
    com_moving_list: List[np.array] = []
    com_reference_list: List[np.array] = []

    # indicies to positions
    for point in points:
        # Reference frame position
        position = interpolated_grid.point_to_position(
            interpolated_grid.get_point(point[0], point[1], point[2]))

        # Tranform to moving frame
        position_moving = transform.apply(position)

        # Interpolate moving map
        interpolated_map_value = moving_xmap.xmap.interpolate_value(
            position_moving)

        # Set original point
        interpolated_grid.set_value(
            point[0], point[1], point[2], interpolated_map_value)

    # interpolated_grid.symmetrize_max()
    interpolated_array = np.array(interpolated_grid)

    interpolated_grid_neg = gridFromTemplate(moving_xmap)
    interpolated_array_neg = np.array(interpolated_grid_neg, copy=False)
    interpolated_array_neg[:, :, :] = -interpolated_array[:, :, :]
    interpolated_grid_neg.symmetrize_max()

    interpolated_grid_pos = gridFromTemplate(moving_xmap)
    interpolated_array_pos = np.array(interpolated_grid_pos, copy=False)
    interpolated_array_pos[:, :, :] = interpolated_array[:, :, :]
    interpolated_grid_pos.symmetrize_max()

    interpolated_grid_sym = gridFromTemplate(moving_xmap)
    interpolated_array_sym = np.array(interpolated_grid_sym, copy=False)
    interpolated_array_sym[:, :, :] = interpolated_array_pos[:,
                                                             :, :] - interpolated_array_neg[:, :, :]

    return Xmap(interpolated_grid_sym)


def referenceSave(template_map_path, xmap, path_to_save):
    # Open Template map
    print(str(template_map_path))
    print(str(path_to_save))
    ccp4 = gemmi.read_ccp4_map(str(template_map_path))
    ccp4.setup()
    # Replace Template map data
    ccp4.grid = xmap.xmap
    ccp4.setup()
    # Save Template Map...
    ccp4.write_ccp4_map(str(path_to_save))


def gridFromTemplate(template):
    interpolated_grid = gemmi.FloatGrid(
        template.xmap.nu,
        template.xmap.nv,
        template.xmap.nw, )
    interpolated_grid.set_unit_cell(template.xmap.unit_cell)
    interpolated_grid.spacegroup = template.xmap.spacegroup
    return interpolated_grid
