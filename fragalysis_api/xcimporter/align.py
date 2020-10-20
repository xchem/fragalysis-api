import glob
import Bio.PDB as bp
import pymol
from pathlib import Path
from scipy import spatial
from fragalysis_api.xcimporter.conversion_pdb_mol import Ligand
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

warnings.simplefilter('ignore', bpp.PDBConstructionWarning)


class Align:

    def __init__(self, directory, pdb_ref=''):

        self.directory = directory
        self._get_ref = pdb_ref

    @property
    def _get_files(self):
        """
        Extracts a list of paths for all PDBs within the given directory.
        :return: list of .pdb file names in directory
        """
        all_files = set(glob.glob(os.path.join(self.directory, '*')))
        txt_files = set(glob.glob(os.path.join(self.directory, '*.txt')))
        return list(all_files - txt_files)

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
        print(str(Path(os.path.join(self.directory, f'{name}{ext}'))))
        map = Xmap.from_file(file=Path(os.path.join(self.directory, f'{name}{ext}')))
        # Cut Map!

        map.resample(xmap=map, transform=transform)
        print(str(Path(os.path.join(out_dir, f'{name}{ext}'))))
        map.save(path=Path(os.path.join(out_dir, f'{name}{ext}')))

    def align(self, out_dir):
        """
        Aligns all pdbs in the pymol object to the pdb_ref.
        :param out_dir: directory to save aligned pdbs in
        :return: saves the pdbs
        """
        # load ref
        input_files = self._get_files
        basenames = [os.path.splitext(os.path.basename(f))[0] for f in input_files]
        crys = [y for y in [x for x in basenames if 'event' not in x] if 'fofc' not in y]

        reference_pdb = Structure.from_file(file=Path(os.path.join(self.directory, f'{self._get_ref}.pdb')))
        for num, name in enumerate(crys):
            fofcs = [j for j in [i for i in basenames if name in i] if '_fofc' in j]
            fofcs2 = [j for j in [i for i in basenames if name in i] if '_2fofc' in j]
            events = [j for j in [i for i in basenames if name in i] if 'event' in j]
            if not name == self._get_ref:
                # Do an alignment + save
                current_pdb = Structure.from_file(file=Path(os.path.join(self.directory, f'{name}.pdb')))
                current_pdb, transform = current_pdb.align_to(other=reference_pdb)
                current_pdb.structure.write_pdb(os.path.join(out_dir, f'{name}_bound.pdb'))
                # Align Xmaps + save!
                print(name)
                for i in fofcs:
                    self.read_reshape_resave(name=i, out_dir=out_dir, ext='.map', transform=transform)

                for i in fofcs2:
                    self.read_reshape_resave(name=i, out_dir=out_dir, ext='.map', transform=transform)

                for i in events:
                    self.read_reshape_resave(name=i, out_dir=out_dir, ext='.ccp4', transform=transform)

            else:
                shutil.copyfile(os.path.join(self.directory, f'{name}.pdb'),
                                os.path.join(out_dir, f'{name}_bound.pdb'))
                shutil.copyfile(os.path.join(self.directory, f'{name}_smiles.txt'),
                                os.path.join(out_dir, f'{name}_smiles.txt'))
                for i in fofcs:
                    shutil.copyfile(os.path.join(self.directory, f'{i}.map'),
                                    os.path.join(out_dir, f'{i}.map'))
                for i in fofcs2:
                    shutil.copyfile(os.path.join(self.directory, f'{i}.map'),
                                    os.path.join(out_dir, f'{i}.map'))
                for i in events:
                    shutil.copyfile(os.path.join(self.directory, f'{i}.ccp4'),
                                    os.path.join(out_dir, f'{i}.ccp4'))


class CutMaps:

    def __init__(self, in_dir, out_dir, monomerize):

        self.in_dir = in_dir
        self.out_dir = out_dir
        self.monomerize = monomerize

    @property
    def _get_files(self):
        """
        Extracts a list of paths for all PDBs within the given directory.
        :return: list of .pdb file names in directory
        """
        all_files = set(glob.glob(os.path.join(self.in_dir, '*')))
        txt_files = set(glob.glob(os.path.join(self.in_dir, '*.txt')))
        return list(all_files - txt_files)

    def cut_maps(self):
        input_files = self._get_files
        basenames = [os.path.splitext(os.path.basename(f))[0] for f in input_files]
        crys = [y for y in [x for x in basenames if 'event' not in x] if 'fofc' not in y]

        for name in crys:
            print(f'Cutting {name}...')
            basepdb = os.path.join(self.in_dir, f'{name}.pdb')
            name = basepdb.replace('_bound', '')
            new = Ligand(name, basepdb, self.out_dir)
            new.hets_and_cons()
            new.remove_nonligands()
            new.find_ligand_names_new()
            fofcmap = os.path.join(self.in_dir, f'{name}_fofc.map')
            fofc2map = os.path.join(self.in_dir, f'{name}_2fofc.map')
            events = [i for i in basenames if f'{name}_event' in i]
            for i, lig_name in enumerate(new.wanted_ligs):
                print(lig_name)
                xyzin, base = new.create_pdb_for_ligand(lig_name, count=i, monomerize=self.monomerize, smiles_file=None, out_dir=self.out_dir, ret2=True)
                print(xyzin)
                fofcout = os.path.join(self.out_dir, f'{base}_fofc.map')
                fofc2out = os.path.join(self.out_dir, f'{base}_2fofc.map')
                # Now cut the maps and copy the files
                cmd = (f"module load ccp4 && mapmask mapin {fofcmap} mapout {fofcout} xyzin {xyzin} << eof\n border 10\n end\n eof")
                os.system(cmd)
                cmd = (f"module load ccp4 && mapmask mapin {fofc2map} mapout {fofc2out} xyzin {xyzin} << eof\n border 10\n end\n eof")
                os.system(cmd)
                for num, j in enumerate(events):
                    eventmap = os.path.join(self.in_dir, f'{j}.ccp4')
                    eventout = os.path.join(self.out_dir, f'{base}_event_{num}.ccp4')
                    cmd = (f"module load ccp4 && mapmask mapin {eventmap} mapout {eventout} xyzin {xyzin} << eof\n border 10\n end\n eof")
                    os.system(cmd)
                # clean-up
                os.remove(xyzin)
                # copy txt file?
                shutil.copyfile(os.path.join(self.in_dir, f'{name}_smiles.txt'),
                                os.path.join(self.out_dir, f'{name}_smiles.txt'))
                shutil.copyfile(os.path.join(self.in_dir, f'{name}.pdb'),
                                os.path.join(self.out_dir, f'{name}.pdb'))


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

    def align_to(self, other):
        # Warning: inplace!
        # Aligns structures usings carbon alphas and transform self into the frame of the other

        ca_self = []
        ca_other = []

        # Get CAs
        for model in self.structure:
            for chain in model:
                for res_self in chain.get_polymer():
                    if 'LIG' in str(res_self):
                        continue

                    try:
                        current_res_id = ResidueID.from_residue_chain(model, chain, res_self)
                        res_other = other.structure[current_res_id.model][current_res_id.chain][current_res_id.insertion][0]
                        print(f'{self.structure}|{res_self}')
                        print(f'{other.structure}|{res_other}')
                        self_ca_pos = res_self["CA"][0].pos
                        other_ca_pos = res_other["CA"][0].pos

                    except:
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
            print(str(atom))
            atom.pos = transform.apply_inverse(atom.pos)

        return self, transform


@dataclasses.dataclass()
class Xmap:
    xmap: gemmi.FloatGrid

    @staticmethod
    def from_file(file):
        ccp4 = gemmi.read_ccp4_map(str(file))
        ccp4.setup()
        return Xmap(ccp4.grid)

    def resample(
            self,
            xmap,
            transform,  # tranfrom FROM the frame of xmap TO the frame of self
            sample_rate: float = 3.0,
    ):
        print('Map XForm!!!')
        unaligned_xmap: gemmi.FloatGrid = self.xmap

        unaligned_xmap_array = numpy.array(unaligned_xmap, copy=False)
        std = numpy.std(unaligned_xmap_array)

        unaligned_xmap_array[:, :, :] = unaligned_xmap_array[:, :, :] / std

        interpolated_values_tuple = ([], [], [], [])

        alignment_positions: typing.Dict[typing.Tuple[int], gemmi.Position] = {
            point: unaligned_xmap.point_to_position(unaligned_xmap.get_point(*point))
            for point, value
            in numpy.ndenumerate(unaligned_xmap_array)
        }

        transformed_positions: typing.Dict[typing.Tuple[int], gemmi.Position] = {
            point: transform.apply(position) for point, position in alignment_positions.items()
        }

        transformed_positions_fractional: typing.Dict[typing.Tuple[int], gemmi.Fractional] = {
            point: unaligned_xmap.unit_cell.fractionalize(pos) for point, pos in transformed_positions.items()}

        interpolated_values: typing.Dict[typing.Tuple[int],
                                         float] = Xmap.interpolate_grid(unaligned_xmap,
                                                                        transformed_positions_fractional)

        interpolated_values_tuple = (interpolated_values_tuple[0] + [index[0] for index in interpolated_values],
                                     interpolated_values_tuple[1] + [index[1] for index in interpolated_values],
                                     interpolated_values_tuple[2] + [index[2] for index in interpolated_values],
                                     interpolated_values_tuple[3] + [interpolated_values[index] for index in
                                                                     interpolated_values],
                                     )

        # Copy data into new grid
        new_grid = xmap.new_grid()
        grid_array = numpy.array(new_grid, copy=False)
        grid_array[interpolated_values_tuple[0:3]] = interpolated_values_tuple[3]
        return Xmap(new_grid)

    def new_grid(self):
        spacing = [self.xmap.nu, self.xmap.nv, self.xmap.nw]
        unit_cell = self.xmap.unit_cell
        grid = gemmi.FloatGrid(spacing[0], spacing[1], spacing[2])
        grid.unit_cell = unit_cell
        grid.spacegroup = self.xmap.spacegroup
        return grid

    @staticmethod
    def interpolate_grid(grid: gemmi.FloatGrid,
                         positions: typing.Dict[typing.Tuple[int],
                                                gemmi.Position]) -> typing.Dict[typing.Tuple[int], float]:
        return {coord: grid.interpolate_value(pos) for coord, pos in positions.items()}

    def to_array(self, copy=True):
        return numpy.array(self.xmap, copy=copy)

    def save(self, path: Path, p1: bool = True):
        ccp4 = gemmi.Ccp4Map()
        ccp4.grid = self.xmap
        if p1:
            ccp4.grid.spacegroup = gemmi.find_spacegroup_by_name("P 1")
        else:
            ccp4.grid.symmetrize_max()
        ccp4.update_ccp4_header(2, True)
        ccp4.write_ccp4_map(str(path))


class Monomerize:

    def __init__(self, directory, outdir):
        self.directory = directory
        self.outdir = outdir
        self.non_ligs = json.load(
            open(os.path.join(os.path.dirname(__file__), "non_ligs.json"), "r")
        )

    def get_filelist(self):
        return glob.glob(os.path.join(self.directory, '*.pdb'))

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

    def save_chain(self, lig, f):
        lig_chain = lig[5]
        name = os.path.splitext(os.path.basename(f))[0] + '_' + str(lig_chain)
        filename = os.path.join(self.outdir, f'{name}_mono.pdb')
        pymol.cmd.load(f, name)
        pymol.cmd.save(filename, f'chain {lig_chain}')
        pymol.cmd.reinitialize()

        return filename

    def process_ligs(self, filename):
        test_block = open(filename, 'r').readlines()
        ligs = self.find_ligs(test_block)
        print(ligs)
        outnames = []
        for lig in ligs:
            o = self.save_chain(lig, filename)
            outnames.append(o)
            if os.path.isfile(filename.replace('.pdb', '_smiles.txt')):
                shutil.copy(filename.replace('.pdb', '_smiles.txt'), o.replace('_mono.pdb', '_smiles.txt'))
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
                print(new_pdb)
                handle.write(new_pdb)

    def monomerize_all(self):
        for f in self.get_filelist():
            print(f)
            outnames = self.process_ligs(f)
            print(outnames)
            self.write_bound(f, outnames)
            for o in outnames:
                if os.path.isfile(o):
                    os.remove(o)
