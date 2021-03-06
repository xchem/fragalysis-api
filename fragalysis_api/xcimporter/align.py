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
import json
import shutil
import gemmi  # Oh boy...
import numpy
import argparse

warnings.simplefilter('ignore', bpp.PDBConstructionWarning)


class Align:

    def __init__(self, directory, pdb_ref='', mono=False):
        '''
        :param directory: Directory path contain pdbs to be aligned.
        :param pdb_ref: The String Reference of a pdb that you want to
            optionally align the input pdbs to. This should be the name of the PDB without
            .pbd and should be located inside the specified directory.
        :param mono: Bool, To indicate if chains have been split into seperate pdb files.
            If True, this will change the behaviour of the alignment to solely align the first
            chain of each structure together, and not consider other chains.
        '''

        self.directory = directory
        self._get_ref = pdb_ref
        self.mono = mono

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
                assert (os.path.isfile(os.path.join(self.directory, pdb_ref + ".pdb")) == 1)
                self.__pdb_ref = pdb_ref
            except AssertionError:
                print('pdb desired as reference does not exists. Default pdb chosen.')
                self.__pdb_ref = self.__best_length_and_resolution([i for i in self._get_files if 'pdb' in i])
        else:
            self.__pdb_ref = self.__best_length_and_resolution([i for i in self._get_files if 'pdb' in i])

    def __best_length_and_resolution(self, pdb_files):
        """
        Find the longest pdb structure with the lowest resolution from all imported files.
        This will be used as the reference pdb for alignment against.
        :param pdb_files: List of pdb filepaths
        :return: str of filename with best .pdb file
        """
        a_df = pd.DataFrame(pdb_files, columns=['file'])
        a_df[['resolution', 'p_len', 'pdb']] = a_df.file.apply(self.__get_length_and_resolution)

        return a_df.sort_values(by=['p_len', 'resolution'], ascending=[False, True]).pdb[0]

    def __get_length_and_resolution(self, file):
        """
        Determine resolution, sequence and length of .pdb file.
        :param file: pdb file path.
        :return: pandas series with resolution, sequence and length and .pdb filename.
        """
        parser = bp.PDBParser()
        ppb = bp.PPBuilder()
        structure = parser.get_structure(os.path.splitext(os.path.basename(file))[0], file)

        seq_len = 0

        for pp in ppb.build_peptides(structure):  # Retrieve length by looping through each chain in the protein
            seq_len += len(pp.get_sequence())

        # using a functions from PDBParser parser class to get the resolution and protein id from the pdb file
        return pd.Series([structure.header['resolution'], seq_len, structure.id])

    def align_to_reference(self, in_file, reference, out_dir):
        '''
        Aligns a single pdb file to a reference and adds it to the specified output directory.
        :param in_file: filepath to corresponding pdb file to align. Accompanying map files should be located within
            the same directory as input file.
        :param reference: the reference pdb file to align to.
        :param out_dir: The desired output directory for the aligned pdb file.
        :return: an aligned pdb file in the output directory with the same name as input.
        '''
        input_files = in_file
        map_list = self._get_maplist
        base_names = os.path.splitext(os.path.basename(input_files))[0]
        crystals = [y for y in [x for x in [base_names] if 'event' not in x] if 'fofc' not in y]
        ref = reference
        dir = self.directory
        mono = self.mono
        for num, name in enumerate(crystals):
            reference_pdb = Structure.from_file(file=Path(ref))
            all_maps = [j for j in map_list if name in j]
            current_pdb = Structure.from_file(file=Path(in_file))
            try:
                current_pdb, transform = current_pdb.align_to(other=reference_pdb, monomerized=mono)
            except Exception as e:
                # Poorly Documented. Use better stuff...
                print(f'{e}')
                continue

            current_pdb.structure.write_pdb(os.path.join(out_dir, f'{name}_bound.pdb'))
            # Align Xmaps + save!
            s = time.time()
            print(all_maps)
            for i in all_maps:
                base, ext = os.path.splitext(os.path.basename(i))
                print(i)
                s2 = time.time()
                map = Xmap.from_file(file=Path(os.path.join(dir, f'{base}{ext}')))
                map.resample(xmap=map, transform=transform)
                map.save(path=Path(os.path.join(out_dir, f'{base}{ext}')))
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
        shutil.copyfile(fn, output)

    def align(self, out_dir):
        """
        Aligns all pdb structures and map files (if any) using gemmi, and save them to a new directory.
        :param out_dir: directory to save aligned pdbs in
        :return: saves the pdbs + transforms map files (if any!)
        """
        # load ref
        input_files = self._get_files
        map_list = self._get_maplist
        crystals = [os.path.splitext(os.path.basename(f))[0] for f in input_files]
        ref = self._get_ref
        dir = self.directory
        mono = self.mono

        s = time.time()
        for num, name in enumerate(crystals):
            reference_pdb = Structure.from_file(file=Path(os.path.join(dir, f'{ref}.pdb')))
            all_maps = [j for j in map_list if name in j]
            if not name == ref:
                # Do an alignment + save
                current_pdb = Structure.from_file(file=Path(os.path.join(dir, f'{name}.pdb')))
                try:
                    current_pdb, transform = current_pdb.align_to(other=reference_pdb, monomerized=mono)
                except Exception as e:
                    # Poorly Documented. Use better stuff...
                    print(f'{e}')
                    continue

                current_pdb.structure.write_pdb(os.path.join(out_dir, f'{name}_bound.pdb'))
                # Align Xmaps + save!
                for i in all_maps:
                    base, ext = os.path.splitext(os.path.basename(i))
                    print(i)
                    # self.read_reshape_resave(name=base, out_dir=out_dir, ext=ext, transform=transform)
                    s2 = time.time()
                    map = Xmap.from_file(file=Path(os.path.join(dir, f'{base}{ext}')))
                    map.resample(xmap=map, transform=transform)
                    map.save(path=Path(os.path.join(out_dir, f'{base}{ext}')))
                    e2 = time.time()
                    print(f'{int(e2 - s2)} seconds to transform map...')

                    e = time.time()
                    print(f'Total Running time: {int(e - s) / 60} minutes.')

            else:
                shutil.copyfile(os.path.join(self.directory, f'{name}.pdb'),
                                os.path.join(out_dir, f'{name}_bound.pdb'))
                if os.path.exists(os.path.join(self.directory, f'{name}_smiles.txt')):
                    shutil.copyfile(os.path.join(self.directory, f'{name}_smiles.txt'),
                                    os.path.join(out_dir, f'{name}_smiles.txt'))
                for i in all_maps:
                    base, ext = os.path.splitext(os.path.basename(i))
                    shutil.copyfile(i, os.path.join(out_dir, f"{base}{ext}"))


# Conor's stuff
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
        rfree = structure.structure.make_mmcif_document()[0].find_loop("_refine.ls_R_factor_R_free")[0]
        return RFree(float(rfree))

    def to_float(self):
        return self.rfree


@dataclasses.dataclass()
class Transform:
    transform: gemmi.Transform
    com_reference: numpy.array
    com_moving: numpy.array

    def apply(self, position: gemmi.Position) -> gemmi.Position:
        rotation_frame_position = gemmi.Position(position[0] - self.com_reference[0],
                                                 position[1] - self.com_reference[1],
                                                 position[2] - self.com_reference[2])
        transformed_vector = self.transform.apply(rotation_frame_position)
        transformed_position = gemmi.Position(transformed_vector[0] + self.com_moving[0],
                                              transformed_vector[1] + self.com_reference[1],
                                              transformed_vector[2] + self.com_reference[2])
        return transformed_position

    def apply_inverse(self, position: gemmi.Position) -> gemmi.Position:
        rotation_frame_position = gemmi.Position(position[0] - self.com_moving[0],
                                                 position[1] - self.com_moving[1],
                                                 position[2] - self.com_moving[2])
        transformed_vector = self.transform.inverse().apply(rotation_frame_position)
        transformed_position = gemmi.Position(transformed_vector[0] + self.com_reference[0],
                                              transformed_vector[1] + self.com_reference[1],
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

    def align_to(self, other, monomerized=False):
        # Warning: inplace!
        # Aligns structures usings carbon alphas and transform self into the frame of the other

        ca_self = []
        ca_other = []

        # Get CAs
        for model in self.structure:
            for chain in model:
                for res_self in chain.get_polymer():
                    if 'LIG' in str(res_self):
                        print('Skipping Ligand...')
                        continue

                    try:
                        current_res_id = ResidueID.from_residue_chain(model, chain, res_self)
                        if monomerized:
                            # print(other.structure[current_res_id.model])
                            # print(len(other.structure[current_res_id.model]))
                            res_other = other.structure[current_res_id.model][0][current_res_id.insertion][0]
                        else:
                            res_other = \
                                other.structure[current_res_id.model][current_res_id.chain][current_res_id.insertion][0]
                        # print(f'{self.structure}|{res_self}')
                        # print(f'{other.structure}|{res_other}')
                        self_ca_pos = res_self["CA"][0].pos
                        other_ca_pos = res_other["CA"][0].pos

                    except:
                        print('Skipping, Residue not found in chain')
                        continue

                    ca_list_self = Transform.pos_to_list(self_ca_pos)
                    ca_list_other = Transform.pos_to_list(other_ca_pos)

                    ca_self.append(ca_list_self)
                    ca_other.append(ca_list_other)

        # Make coord matricies
        matrix_self = numpy.array(ca_self)
        matrix_other = numpy.array(ca_other)

        # Find means
        mean_self = numpy.mean(matrix_self, axis=0)
        mean_other = numpy.mean(matrix_other, axis=0)

        # demaen
        de_meaned_self = matrix_self - mean_self
        de_meaned_other = matrix_other - mean_other

        # Align
        rotation, rmsd = spatial.transform.Rotation.align_vectors(de_meaned_self, de_meaned_other)

        # Get transform
        vec = numpy.array([0.0, 0.0, 0.0])
        # Transform is from other frame to self frame
        transform = Transform.from_translation_rotation(vec,
                                                        rotation,
                                                        mean_other,
                                                        mean_self,
                                                        )

        # Transform positions
        for atom in self.all_atoms():
            # print(str(atom))
            atom.pos = transform.apply_inverse(atom.pos)

        return self, transform


class Monomerize:

    def __init__(self, directory, outdir):
        '''
        :param directory: Path to folder containing .pdb files to monomerize
        :param outdir: Output folder where monomerized pdbs will be saved
        :return: Initialises the Monomerize class that has the ability to split pdb files by chain.
        '''
        self.directory = directory
        self.outdir = outdir
        self.non_ligs = json.load(
            open(os.path.join(os.path.dirname(__file__), "non_ligs.json"), "r")
        )

    def get_filelist(self):
        '''
        Get a list of pdb files within the input directory
        :return: Returns a list of pdb files
        '''
        return glob.glob(os.path.join(self.directory, '*.pdb'))

    def get_maplist(self):
        '''
        Get a list of .map or .ccp4 files from within the input directory
        :return: Returns a list of filepaths corresponding to said volume density files.
        '''
        map_files = glob.glob(os.path.join(self.directory, '*.map'))
        cpp4_files = glob.glob(os.path.join(self.directory, '*.ccp4'))
        return map_files + cpp4_files

    def split_chains(self, f):
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
                print(chain_dists)
                temp_structure[0][j].name = min(chain_dists, key=chain_dists.get)

            leftover_chains = list(set(nonHOH_chains) - set(i))
            # Remove remaining chains
            for j in leftover_chains:
                if j == i:
                    # Just in case.
                    continue
                else:
                    temp_structure[0].remove_chain(j)
                    print([x.name for x in temp_structure[0]])

            # Rename Chain to corresponding chain then save!
            name = os.path.splitext(os.path.basename(f))[0] + '_' + str(i)
            filename = os.path.join(self.outdir, f'{name}_mono.pdb')
            print(f'Writing to: {filename}')
            temp_structure.write_pdb(filename)
            filenames.append(filename)

        return filenames

    def process_pdb(self, filename, maplist):
        '''
        Processes a pdb file.
        Firstly it splits the pdb file according to chains and then makes a copy all accompanying files (e.g. map files) to share the same name as the new split pdb file.
        :param filename: The filepath of the pdb file to process
        :param maplist: A list of mapfiles to be copied for the number of chains the filename may have.
        :return: A list of filenames that were generated from the input pdb.
        '''
        out_names = self.split_chains(filename)
        for o in out_names:
            if os.path.isfile(filename.replace('.pdb', '_smiles.txt')):
                shutil.copy(filename.replace('.pdb', '_smiles.txt'), o.replace('_mono.pdb', '_smiles.txt'))
            base = os.path.splitext(os.path.basename(filename))[0]
            new = os.path.splitext(os.path.basename(o))[0].replace('_mono', '')
            allmaps = [j for j in maplist if base in j]
            for map in allmaps:
                if os.path.isfile(map):
                    mapbase = os.path.basename(map)
                    shutil.copy(map, os.path.join(self.outdir, mapbase.replace(base, new)))
        return out_names

    def write_bound(self, inname, outname):
        '''
        Appends header to pdb file and writes to new location
        :param inname: input pdb file
        :param outname: filepath of desired output
        :return: Formally returns nothing, but a file will be written to outname
        '''
        with open(inname, 'r') as handle:
            switch = 0
            header_front, header_end = [], []

            for line in handle:

                if line.startswith('ATOM'): switch = 1

                if line.startswith('HETATM'): switch = 2

                if switch == 0 and not line.startswith('REMARK 350'): header_front.append(line)

                if (switch == 2) and not line.startswith('HETATM'): header_end.append(line)

        for o in outname:
            newfile_contents = open(o, 'r').readlines()

            with open(o.replace('_mono.pdb', '.pdb'), 'w') as handle:
                remark = ['REMARK warning: chains may be ommitted for alignment\n']
                new_pdb = ''.join(remark + header_front + newfile_contents)
                # print(new_pdb)
                handle.write(new_pdb)

    def monomerize_single(self, file):
        '''
        Convert a pdb file into it's constituent chains
        :param file: Input pdb filepath
        :return: Formally nothing, but will create n number of monomerized pdbs per number of AA containing chains.
        '''
        maplist = self.get_maplist()
        outnames = self.process_pdb(file, maplist)
        print(outnames)
        self.write_bound(file, outnames)
        for o in outnames:
            if os.path.isfile(o):
                os.remove(o)

    def monomerize_all(self):
        '''
        Converts a set of pdbs into constituent chains
        :return: Formally nothing, but will separate pdb files according to chains.
        '''
        for f in self.get_filelist():
            self.monomerize_single(file=f)


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
