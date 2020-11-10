import glob
import time

import Bio.PDB as bp
import pymol
from pathlib import Path
from scipy import spatial
# from fragalysis_api import Ligand
from pandda_gemmi.pandda_types import *
import typing
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

        self.directory = directory
        self._get_ref = pdb_ref
        self.mono = mono

    @property
    def _get_files(self):
        """
        Extracts a list of paths for all PDBs within the given directory.
        :return: list of .pdb file names in directory
        """
        all_files = set(glob.glob(os.path.join(self.directory, '*')))
        txt_files = set(glob.glob(os.path.join(self.directory, '*.txt')))
        return list(all_files - txt_files)

    @property
    def _get_maplist(self):
        all_files = set(glob.glob(os.path.join(self.directory, '*')))
        txt_files = set(glob.glob(os.path.join(self.directory, '*.txt')))
        pdb_files = set(glob.glob(os.path.join(self.directory, '*.pdb')))
        return list(all_files - txt_files - pdb_files)

    def _load_objs(self):
        """
        Loads each pdb into the PyMol instance/object.
        :type: object
        :return: PyMol object with .pdb protein structure file loaded
        """
        # Looping through each pdb file in the directory and loading them into the cmd
        for num, file in enumerate(self._get_files):
            pymol.cmd.load(file, os.path.splitext(os.path.basename(file))[0])

        # deal with files that have no conect records
        pymol.cmd.set('pdb_conect_all', 'on')

        return pymol.cmd

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
        :return PyMol instance with reference object assigned as reference property:
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
        :param pdb_files:
        :return: str of filename with best .pdb file
        """
        a_df = pd.DataFrame(pdb_files, columns=['file'])
        a_df[['resolution', 'p_len', 'pdb']] = a_df.file.apply(self.__get_length_and_resolution)

        return a_df.sort_values(by=['p_len', 'resolution'], ascending=[False, True]).pdb[0]

    def __get_length_and_resolution(self, file):
        """
        Determine resolution, sequence and length of .pdb file.
        :param file:
        :return pandas series with resolution, sequence and length and .pdb filename:
        """
        parser = bp.PDBParser()
        ppb = bp.PPBuilder()
        structure = parser.get_structure(os.path.splitext(os.path.basename(file))[0], file)

        seq_len = 0

        for pp in ppb.build_peptides(structure):  # Retrieve length by looping through each chain in the protein
            seq_len += len(pp.get_sequence())

        # using a functions from PDBParser parser class to get the resolution and protein id from the pdb file
        return pd.Series([structure.header['resolution'], seq_len, structure.id])

    def __get_pdb_file(self, pdb_name):
        """
        Returns the pdb file path for a given pdb name.
        :param pdb_name: name of pdb
        :return: file path to the pdb with the given name
        """

        for file in self._get_files:
            if pdb_name == os.path.splitext(os.path.basename(file))[0]:
                return file

    def __get_header(self, pdb_file):
        """
        Identifies the section of a PDB which contains the headers ATOM/HETATM
        :param pdb_file: The pdb to acquire the header locations of
        :return: front locations of the ATOM/HETATM headers in the given pdb, end locations of the ATOM/HETATM headers in the given pdb       
        """
        with open(pdb_file) as handle:
            switch = 0
            header_front, header_end = [], []

            for line in handle:

                if line.startswith('ATOM'): switch = 1

                if line.startswith('HETATM'): switch = 2

                if switch == 0: header_front.append(line)

                if (switch == 2) and not line.startswith('HETATM'): header_end.append(line)

        return header_front, header_end

    def _save_align(self, name, pdb_base, out_dir):
        """
        Saves modified pdb as pdb file again. It also ensures the pdb header is kept in the new file.
        :param name: name of pdb
        :param pdb_base: coordinates of atom in the pdb format
        :param out_dir: directory to save new pdb file in
        :return: a saved pdb file
        """
        if not os.path.exists(out_dir):  # Creating output directory if it doesn't already exist
            os.makedirs(out_dir)

        pdb_file = self.__get_pdb_file(name)
        header_front = self.__get_header(pdb_file)[0]

        with open(os.path.join(out_dir, f'{name}_bound.pdb'), 'w') as handle:

            new_pdb = header_front + [pdb_base]

            for line in new_pdb:
                handle.write(line)

    def read_reshape_resave(self, name, out_dir, ext, transform):
        map = Xmap.from_file(file=Path(os.path.join(self.directory, f'{name}{ext}')))
        # Cut Map!
        s2 = time.time()
        map.resample(xmap=map, transform=transform)
        e2 = time.time()
        print(f'{int(e2 - s2) / 60} minutes ({int(e2 - s2)} seconds) taken to transform map...')
        s = time.time()
        map.save(path=Path(os.path.join(out_dir, f'{name}{ext}')))
        e = time.time()
        print(f'{int(e - s) / 60} minutes ({int(e - s)} seconds) taken to save map...')

    def align_to_reference(self, in_file, reference, out_dir):
        '''

        Parameters
        ----------
        in_file: path to input pdb files (str)
        reference: path to input pdb file (str)
        out_dir: desired output directory (str)

        Returns
        -------

        '''
        input_files = in_file
        map_list = self._get_maplist
        base_names = [os.path.splitext(os.path.basename(f))[0] for f in input_files]
        crystals = [y for y in [x for x in base_names if 'event' not in x] if 'fofc' not in y]
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


    def align(self, out_dir):
        """
        Aligns all pdbs in the pymol object to the pdb_ref.
        :param monomerized: Logical, whether chains have been split into single PDBs
        :param out_dir: directory to save aligned pdbs in
        :return: saves the pdbs + transforms map files (if any!)
        """
        # load ref
        input_files = self._get_files
        map_list = self._get_maplist
        base_names = [os.path.splitext(os.path.basename(f))[0] for f in input_files]
        crystals = [y for y in [x for x in base_names if 'event' not in x] if 'fofc' not in y]
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
                    #self.read_reshape_resave(name=base, out_dir=out_dir, ext=ext, transform=transform)
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
                            #print(other.structure[current_res_id.model])
                            #print(len(other.structure[current_res_id.model]))
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
        self.directory = directory
        self.outdir = outdir
        self.non_ligs = json.load(
            open(os.path.join(os.path.dirname(__file__), "non_ligs.json"), "r")
        )

    def get_filelist(self):
        return glob.glob(os.path.join(self.directory, '*.pdb'))

    def get_maplist(self):
        all_files = set(glob.glob(os.path.join(self.directory, '*')))
        txt_files = set(glob.glob(os.path.join(self.directory, '*.txt')))
        pdb_files = set(glob.glob(os.path.join(self.directory, '*.pdb')))
        return list(all_files - txt_files - pdb_files)

    def find_ligs(self, pdb_lines):
        """
        Finds list of ligands contained in the structure, including
        """
        all_ligands = []  # all ligands go in here, including solvents and ions
        wanted_ligs = []
        for line in pdb_lines:
            if line.startswith("HETATM"):
                all_ligands.append(line)

        for lig in all_ligands:
            if (
                    lig.split()[3][-3:] not in self.non_ligs
            ):  # this takes out the solvents and ions a.k.a non-ligands
                wanted_ligs.append(lig[16:20] + lig[20:26])
                # print(lig[16:20].strip() + lig[20:26])

        wanted_ligs = list(set(wanted_ligs))

        return wanted_ligs

    def split_chains(self, f):
        filenames = []
        base_structure = gemmi.read_structure(f)
        base_models = base_structure[0]
        all_chains = [x.name for x in base_models if not x.name == 'S']
        chain_centers = {}
        for i in all_chains:
            chain_centers[i] = get_chain_center(chain_name=i, file=f)

        print(chain_centers)

        chain_names = [x.name for x in base_models if x.calculate_mass() > 5000 and not x.name == 'S']
        print(chain_names)
        alt_chains = [x.name for x in base_models if x.calculate_mass() <= 5000 and not x.name == 'S']
        print(alt_chains)

        for i in chain_names:
            # For each chain, convert all ligands,
            temp_structure = gemmi.read_structure(f)
            for j in alt_chains:
                chain_dists = {}
                for z in chain_names:
                    chain_dists[z] = chain_centers[j].dist(chain_centers[z])
                print(chain_dists)
                temp_structure[0][j].name = min(chain_dists, key=chain_dists.get)

            # Flatten to single chain
            print([x.name for x in temp_structure[0]])
            # temp_structure.merge_chain_parts()
            print([x.name for x in temp_structure[0]])
            leftover_chains = [x.name for x in temp_structure[0] if not x.name == i]
            print(leftover_chains)

            # Remove remaining chains
            for j in leftover_chains:
                if j == i:
                    continue
                else:
                    temp_structure[0].remove_chain(j)
                    print([x.name for x in temp_structure[0]])

            # Rename Chain to corresponding chain then save!
            name = os.path.splitext(os.path.basename(f))[0] + '_' + str(i)
            filename = os.path.join(self.outdir, f'{name}_mono.pdb')
            print(f'Writing to: {filename}')
            temp_structure.write_pdb(filename)
            #temp_structure.write_minimal_pdb(filename)
            filenames.append(filename)

        return filenames

    def process_pdb(self, filename, maplist):
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

    def save_chain(self, lig, f):

        # This is wrong...
        # For each ligand, set and assign the closest chain to it.
        # Then after all the ligands are set correct
        # Save each individual chain!
        # Figure out what lig and f are...
        print(lig)
        lig_chain = lig[5]

        base_structure = gemmi.read_structure(f)
        base_models = base_structure[0]

        ligchain = base_models.find_chain(lig_chain)
        ligchain_pos = ligchain[0][0].pos

        chain_dists = {}
        for chain in base_models:
            chain_name = chain.name

            # Defo the wrong way to go about this!!
            if chain.has_subchains_assigned() and chain.calculate_mass() > 5000:
                chain_temp = gemmi.read_structure(f)[0]

                for k in [j.name for j in chain_temp]:
                    if not k == chain_name:
                        chain_temp.remove_chain(k)

                chain_dists[chain_name] = ligchain_pos.dist(chain_temp.calculate_center_of_mass())

        if len(chain_dists) > 0: # I don't think this is sufficient...
            # Remove remaining chain parts
            closest_chain = min(chain_dists, key=chain_dists.get)
            ligchain.name = closest_chain
            base_structure.merge_chain_parts()
            for x in [j.name for j in base_models if not j.name == closest_chain]:
                base_models.remove_chain(x)

        # Rename Chain to corresponding chain then save!
        name = os.path.splitext(os.path.basename(f))[0] + '_' + str(lig_chain)
        filename = os.path.join(self.outdir, f'{name}_mono.pdb')
        print(f'Writing to: {filename}')
        base_structure.write_pdb(filename)

        return filename

    def process_ligs(self, filename, maplist):
        test_block = open(filename, 'r').readlines()
        ligs = self.find_ligs(test_block)
        print(ligs)
        outnames = []
        for lig in ligs:
            o = self.save_chain(lig, filename)
            outnames.append(o)
            if os.path.isfile(filename.replace('.pdb', '_smiles.txt')):
                shutil.copy(filename.replace('.pdb', '_smiles.txt'), o.replace('_mono.pdb', '_smiles.txt'))

            # Find Maps...
            base = os.path.splitext(os.path.basename(filename))[0]
            new = os.path.splitext(os.path.basename(o))[0].replace('_mono', '')
            allmaps = [j for j in maplist if base in j]
            for map in allmaps:
                if os.path.isfile(map):
                    mapbase = os.path.basename(map)
                    shutil.copy(map, os.path.join(self.outdir, mapbase.replace(base, new)))

        return outnames

    def write_bound(self, inname, outname):
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
                #print(new_pdb)
                handle.write(new_pdb)

    def monomerize_single(self, file):
        maplist = self.get_maplist()
        outnames = self.process_pdb(file, maplist)
        self.write_bound(file, outnames)
        for o in outnames:
            if os.path.isfile(o):
                os.remove(o)

    def monomerize_all(self):
        maplist = self.get_maplist()
        for f in self.get_filelist():
            # print(f)
            #outnames = self.process_ligs(f, maplist)
            outnames = self.process_pdb(f, maplist)
            print(outnames)
            self.write_bound(f, outnames)
            for o in outnames:
                if os.path.isfile(o):
                    os.remove(o)


def get_chain_center(chain_name, file):
    structure = gemmi.read_pdb(file)[0]
    names = [x.name for x in structure]
    for i in names:
        if i == chain_name:
            continue
        else:
            structure.remove_chain(i)

    return structure.calculate_center_of_mass()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m", "--monomerize", action="store_true", default=False, help="Monomerize input"
    )
    parser.add_argument(
        "-a", "--align", action="store_true", default=False, help="Align input"
    )
    parser.add_argument("-t", "--target", help="Target name", required=True)
    args = vars(parser.parse_args())
    monomerize = args["monomerize"]
    align = args['align']
    tar = args['target']
    m = Monomerize(f'/dls/science/groups/i04-1/fragprep/input_test/{tar}/', '/dls/science/groups/i04-1/software/tyler/monotest')
    m.monomerize_all()
    if align:
        a = Align('/dls/science/groups/i04-1/software/tyler/monotest', mono=monomerize)
        a.align('/dls/science/groups/i04-1/software/tyler/tmptest')
