import shutil
import unittest
import os
from pathlib import Path

from fragalysis_api import Align, set_up
from glob import glob
from shutil import rmtree

from fragalysis_api.xcimporter.align import Structure


class AlignTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dir_input = os.path.join('tests', 'data_for_tests')
        cls.dir_output = os.path.join('tests', 'data_for_tests')


class EasyAlign(AlignTest):

    @classmethod
    def setUpClass(cls):
        super(EasyAlign, cls).setUpClass()
        cls.align_obj = Align(os.path.join(
            cls.dir_input, 'examples_to_test0'), pdb_ref='')
        cls.align_obj_w_ref = Align(os.path.join(
            cls.dir_input, 'examples_to_test0'), pdb_ref='6hi3')
        cls.align_obj_w_wrong_ref = Align(os.path.join(
            cls.dir_input, 'examples_to_test0'), pdb_ref='wrong_pdb')
        cls.align_obj_w_maps = Align(os.path.join(
            cls.dir_input, 'examples_to_test5'), pdb_ref='', rrf=False)
        cls.align_obj_w_maps_rrf = Align(os.path.join(
            cls.dir_input, 'examples_to_test5'), pdb_ref='', rrf=True)

    @classmethod
    def tearDownClass(cls):
        # [rmtree(a_path) for a_path in cls.test_path_list]
        pass

    def test_get_files(self):
        files = ['6epu.pdb', '6epv.pdb', '6epx.pdb', '6hi3.pdb']
        whole_files = sorted(
            [os.path.join(self.align_obj.directory, x) for x in files])
        self.assertCountEqual(self.align_obj._get_files, whole_files)

    def test_get_mapfiles(self):
        self.assertCountEqual(self.align_obj._get_maplist, [])

    def test_a_align_with_maps(self):
        dir = os.path.join('tests', 'data_for_tests', 'tmp_map')
        if not os.path.exists(dir):
            os.makedirs(dir)
        self.align_obj_w_maps.align(out_dir=dir)
        map_test_cases = ['Mpro-x0978_bound.pdb', 'Mpro-x0981_bound.pdb', 'Mpro-x2097_bound.pdb',
                          'Mpro-x2119_bound.pdb', 'Mpro-x2097_event.ccp4', 'Mpro-x2119_event.ccp4']
        map_test_exists = [os.path.exists(
            os.path.join(dir, x)) for x in map_test_cases]
        [self.assertTrue(x) for x in map_test_exists]

    def test_b_align_to_reference_with_map(self):
        dir = os.path.join('tests', 'data_for_tests', 'tmp_map_single')
        file = os.path.join('tests', 'data_for_tests',
                            'examples_to_test5', 'Mpro-x2119.pdb')
        ref = os.path.join('tests', 'data_for_tests',
                           'examples_to_test5', 'Mpro-x0978.pdb')
        if not os.path.exists(dir):
            os.makedirs(dir)
        self.align_obj_w_maps.align_to_reference(in_file=file,
                                                 reference_pdb=ref,
                                                 out_dir=dir)
        test_cases = ['Mpro-x2119_bound.pdb', 'Mpro-x2119_event.ccp4']
        test_exists = [os.path.exists(os.path.join(dir, x))
                       for x in test_cases]
        [self.assertTrue(x) for x in test_exists]

    def test_align_rrf(self):
        dir = os.path.join('tests', 'data_for_tests', 'tmp_rrf')
        if not os.path.exists(dir):
            os.makedirs(dir)
        self.align_obj_w_maps_rrf.align(out_dir=dir)
        rrf_test_cases = ['Mpro-x0978_A_bound.pdb', 'Mpro-x0981_A_bound.pdb', 'Mpro-x2097_A_bound.pdb',
                          'Mpro-x2119_A_bound.pdb', 'Mpro-x2097_A_event.ccp4', 'Mpro-x2119_A_event.ccp4']
        rrf_test_exists = [os.path.exists(
            os.path.join(dir, x)) for x in rrf_test_cases]
        [self.assertTrue(x) for x in rrf_test_exists]

    def test_z_structure_class_tests(self):
        path = os.path.join('tests', 'data_for_tests',
                            'examples_to_test5', 'Mpro-x2119.pdb')
        struc = Structure.from_file(Path(path))
        self.assertEqual(struc.rfree().to_float(), 0.2173)
        self.assertEqual(len(struc.residue_ids()), 304)
        self.assertEqual(len([i for i in struc.protein_atoms()]), 2374)
        self.assertEqual(len([i for i in struc.all_atoms()]), 2732)

    def test_g_conversion_pdb_mol(self):
        dir = os.path.join('tests', 'data_for_tests', 'conv')
        dir2 = os.path.join(dir, 'target')
        insmiles = os.path.join('tests', 'data_for_tests',
                                'examples_to_test5', 'Mpro-x2119_smiles.txt')
        # Requires aligned pdb file and text file...
        for mono in [False, True]:
            for covalent in [False, True]:
                for biomol in [None, os.path.join('tests', 'data_for_tests', 'example_header.txt')]:
                    for smiles_file in [None, os.path.abspath(insmiles)]:
                        if not os.path.exists(dir):
                            os.makedirs(dir)
                        if not os.path.exists(dir2):
                            os.makedirs(dir2)

                        if mono:
                            inpdb = os.path.join(
                                'tests', 'data_for_tests', 'tmp_rrf', 'Mpro-x2119_A_bound.pdb')
                            final_path_dir = os.path.join('tests', 'data_for_tests', 'conv', 'target', 'aligned',
                                                          'target-Mpro-x2119_0A')
                            test_cases = ['target-Mpro-x2119_0A.mol',
                                          'target-Mpro-x2119_0A_apo-solv.pdb',
                                          'target-Mpro-x2119_0A_meta.csv',
                                          'target-Mpro-x2119_0A.pdb',
                                          'target-Mpro-x2119_0A_apo.pdb',
                                          'target-Mpro-x2119_0A.sdf',
                                          'target-Mpro-x2119_0A_bound.pdb',
                                          'target-Mpro-x2119_0A_apo-desolv.pdb',
                                          'target-Mpro-x2119_0A_event.ccp4']
                        else:
                            inpdb = os.path.join(
                                'tests', 'data_for_tests', 'tmp_map', 'Mpro-x2119_bound.pdb')
                            final_path_dir = os.path.join('tests', 'data_for_tests', 'conv', 'target', 'aligned',
                                                          'target-Mpro-x2119_0')
                            test_cases = ['target-Mpro-x2119_0.mol',
                                          'target-Mpro-x2119_0_apo-solv.pdb',
                                          'target-Mpro-x2119_0_meta.csv',
                                          'target-Mpro-x2119_0.pdb',
                                          'target-Mpro-x2119_0_apo.pdb',
                                          'target-Mpro-x2119_0.sdf',
                                          'target-Mpro-x2119_0_bound.pdb',
                                          'target-Mpro-x2119_0_apo-desolv.pdb',
                                          'target-Mpro-x2119_0_event.ccp4']

                        set = set_up(
                            target_name='target',
                            infile=os.path.abspath(inpdb),
                            out_dir=dir,
                            rrf=mono,
                            smiles_file=smiles_file,
                            biomol=biomol,
                            covalent=covalent
                        )
                        # Assert the Files were made:
                        test_exists = [os.path.exists(os.path.join(
                            final_path_dir, x)) for x in test_cases]
                        [self.assertTrue(x) for x in test_exists]
                        # Clean-up the files.
                        shutil.rmtree(dir2)

    def test_write_ref(self):
        fp = os.path.join('tests', 'data_for_tests')
        self.align_obj_w_ref.write_align_ref(fp)
        self.assertTrue(os.path.exists(fp))

    def test_get_ref_automatically(self):
        """
        Tests it correctly automatically retrieves the best pdb to use as reference for alignments
        """
        self.assertEqual(self.align_obj._get_ref, '6epv')

    def test_get_ref_when_input_ref(self):
        """
        Tests it correctly assigns the input ref as the reference pdb for alignments
        """
        self.assertEqual(self.align_obj_w_ref._get_ref, '6hi3')

    def test_error_for_inserting_wrong_pdb_ref(self):
        self.assertNotEqual(self.align_obj_w_wrong_ref._get_ref, 'wrong_pdb')
        self.assertEqual(self.align_obj_w_wrong_ref._get_ref, '6epv')

    def test_aligns_saved_correctly(self):
        # align_obj = Align(ATAD2_directory)
        # self.assertIsNotNone(align_obj._save_align)
        pass

    def test_aligns_are_correct_done(self):
        pass


if __name__ == '__main__':
    unittest.main()
