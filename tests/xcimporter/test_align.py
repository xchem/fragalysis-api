import shutil
import unittest
import os
from pathlib import Path

from fragalysis_api import Align, Monomerize, set_up
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
        cls.align_obj = Align(os.path.join(cls.dir_input, 'examples_to_test0'), pdb_ref='')
        cls.align_obj_w_ref = Align(os.path.join(cls.dir_input, 'examples_to_test0'), pdb_ref='6hi3')
        cls.align_obj_w_wrong_ref = Align(os.path.join(cls.dir_input, 'examples_to_test0'), pdb_ref='wrong_pdb')
        cls.align_obj_w_maps = Align(os.path.join(cls.dir_input, 'examples_to_test5'), pdb_ref='')

        # Do a monomerize
        out = os.path.join(cls.dir_output, f'mono')
        if not os.path.isdir(out):
            os.makedirs(out)
        cls.monomerize_obj = Monomerize(directory=os.path.join(cls.dir_input, 'examples_to_test5'), outdir=out)

    @classmethod
    def tearDownClass(cls):
        # [rmtree(a_path) for a_path in cls.test_path_list]
        pass

    def test_get_files(self):
        files = ['6epu.pdb', '6epv.pdb', '6epx.pdb', '6hi3.pdb']
        whole_files = sorted([os.path.join(self.align_obj.directory, x) for x in files])
        self.assertCountEqual(self.align_obj._get_files, whole_files)

    def test_get_mapfiles(self):
        self.assertCountEqual(self.align_obj._get_maplist, [])

    def test_a_monomerize_single(self):
        dir = os.path.join('tests', 'data_for_tests', 'mono')
        if not os.path.exists(dir):
            os.makedirs(dir)
        self.monomerize_obj.monomerize_single(
            file=os.path.join('tests', 'data_for_tests', 'examples_to_test5', 'Mpro-x2097.pdb'))
        single_test_cases = ['Mpro-x2097_A.pdb', 'Mpro-x2097_A_event.ccp4']
        single_test_exists = [os.path.exists(os.path.join('tests', 'data_for_tests', 'mono', x)) for x in
                              single_test_cases]
        print(single_test_exists)
        [self.assertTrue(x) for x in single_test_exists]

    def test_b_monomerize_all(self):
        dir = os.path.join('tests', 'data_for_tests', 'mono')
        if not os.path.exists(dir):
            os.makedirs(dir)
        self.monomerize_obj.monomerize_all()
        all_test_cases = ['Mpro-x0978_A.pdb', 'Mpro-x0981_A.pdb',
                          'Mpro-x2097_A.pdb', 'Mpro-x2097_A_event.ccp4',
                          'Mpro-x2119_A.pdb', 'Mpro-x2119_A_event.ccp4']
        all_test_exists = [os.path.exists(os.path.join('tests', 'data_for_tests', 'mono', x)) for x in all_test_cases]
        print(all_test_exists)
        [self.assertTrue(x) for x in all_test_exists]

    def test_c_align_from_mono(self):
        print(os.path.exists(os.path.join('tests', 'data_for_tests', 'mono')))
        a = Align(os.path.join('tests', 'data_for_tests', 'mono'), pdb_ref='', mono=True)
        dir = os.path.join('tests', 'data_for_tests', 'tmp')
        if not os.path.exists(dir):
            os.makedirs(dir)
        a.align(out_dir=dir)

        a_test_cases = ['Mpro-x0978_A_bound.pdb', 'Mpro-x0981_A_bound.pdb', 'Mpro-x2097_A_bound.pdb',
                        'Mpro-x2119_A_bound.pdb', 'Mpro-x2097_A_event.ccp4', 'Mpro-x2119_A_event.ccp4']
        a_test_exists = [os.path.exists(os.path.join(dir, x)) for x in a_test_cases]
        print(a_test_exists)
        [self.assertTrue(x) for x in a_test_exists]

        dir = os.path.join('tests', 'data_for_tests', 'tmpa2')
        if not os.path.exists(dir):
            os.makedirs(dir)

        a2 = Align(os.path.join('tests', 'data_for_tests', 'mono'), pdb_ref='Mpro-x2119', mono=True)
        a2.align(out_dir=dir)

        a2_test_cases = ['Mpro-x0978_A_bound.pdb', 'Mpro-x0981_A_bound.pdb', 'Mpro-x2097_A_bound.pdb',
                         'Mpro-x2119_A_bound.pdb', 'Mpro-x2097_A_event.ccp4', 'Mpro-x2119_A_event.ccp4']
        a2_test_exists = [os.path.exists(os.path.join(dir, x)) for x in a2_test_cases]
        print(a2_test_exists)
        [self.assertTrue(x) for x in a2_test_exists]

        # Also test if error state ruins things
        # error = Align(os.path.join('tests', 'data_for_tests', 'mono'), pdb_ref='Mpro-x2119', mono=False)

    def test_d_align_with_maps(self):
        dir = os.path.join('tests', 'data_for_tests', 'tmp_map')
        if not os.path.exists(dir):
            os.makedirs(dir)
        self.align_obj_w_maps.align(out_dir=dir)
        map_test_cases = ['Mpro-x0978_bound.pdb', 'Mpro-x0981_bound.pdb', 'Mpro-x2097_bound.pdb',
                          'Mpro-x2119_bound.pdb', 'Mpro-x2097_event.ccp4', 'Mpro-x2119_event.ccp4']
        map_test_exists = [os.path.exists(os.path.join(dir, x)) for x in map_test_cases]
        print(map_test_exists)
        [self.assertTrue(x) for x in map_test_exists]

    def test_e_align_to_reference_with_map(self):
        dir = os.path.join('tests', 'data_for_tests', 'tmp_map_single')
        file = os.path.join('tests', 'data_for_tests', 'examples_to_test5', 'Mpro-x2119.pdb')
        ref = os.path.join('tests', 'data_for_tests', 'examples_to_test5', 'Mpro-x0978.pdb')
        if not os.path.exists(dir):
            os.makedirs(dir)
        self.align_obj_w_maps.align_to_reference(in_file=file,
                                                 reference=ref,
                                                 out_dir=dir)
        test_cases = ['Mpro-x2119_bound.pdb', 'Mpro-x2119_event.ccp4']
        test_exists = [os.path.exists(os.path.join(dir, x)) for x in test_cases]
        [self.assertTrue(x) for x in test_exists]

    def test_z_structure_class_tests(self):
        path = os.path.join('tests', 'data_for_tests', 'examples_to_test5', 'Mpro-x2119.pdb')
        struc = Structure.from_file(Path(path))
        self.assertEqual(struc.rfree().to_float(), 0.2173)
        self.assertEqual(len(struc.residue_ids()), 304)
        self.assertEqual(len([i for i in struc.protein_atoms()]), 2374)
        self.assertEqual(len([i for i in struc.all_atoms()]), 2732)

    def test_f_monomerize_multichain(self):
        dir = os.path.join('tests', 'data_for_tests', 'monomc')
        if not os.path.exists(dir):
            os.makedirs(dir)

        in_dir = os.path.join('tests', 'data_for_tests', 'examples_to_test6')
        mono = Monomerize(directory=in_dir, outdir=dir)
        mono.monomerize_all()
        test_cases = [
            'mArh-x0091_A.pdb',
            'mArh-x0104_A.pdb',
            'mArh-x0128_A.pdb',
            'mArh-x0142_A.pdb',
            'mArh-x0158_A.pdb',
            'mArh-x0091_B.pdb',
            'mArh-x0104_B.pdb',
            'mArh-x0128_B.pdb',
            'mArh-x0142_B.pdb',
            'mArh-x0158_B.pdb'
        ]
        test_exists = [os.path.exists(os.path.join(dir, x)) for x in test_cases]
        [self.assertTrue(x) for x in test_exists]

    def test_g_conversion_pdb_mol(self):
        dir = os.path.join('tests', 'data_for_tests', 'conv')
        dir2 = os.path.join(dir, 'target')
        insmiles = os.path.join('tests', 'data_for_tests', 'examples_to_test5', 'Mpro-x2119_smiles.txt')
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
                            inpdb = os.path.join('tests', 'data_for_tests', 'tmp', 'Mpro-x2119_A_bound.pdb')
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
                            inpdb = os.path.join('tests', 'data_for_tests', 'tmp_map', 'Mpro-x2119_bound.pdb')
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
                            monomerize=mono,
                            smiles_file=smiles_file,
                            biomol=biomol,
                            covalent=covalent
                        )
                        # Assert the Files were made:
                        test_exists = [os.path.exists(os.path.join(final_path_dir, x)) for x in test_cases]
                        print(mono)
                        print(smiles_file)
                        print(biomol)
                        print(covalent)
                        print(test_cases)
                        print(test_exists)
                        print([i for i in os.walk(final_path_dir)])
                        [self.assertTrue(x) for x in test_exists]
                        # Clean-up the files.
                        shutil.rmtree(dir2)


    def test_write_ref(self):
        fp = os.path.join('tests', 'data_for_tests' 'reference.pdb')
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
