import unittest
import os
from fragalysis_api import Align, Monomerize
from glob import glob
from shutil import rmtree


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
        out = os.path.join(cls.dir_output, f'mono/')
        if not os.path.isdir(out):
            os.makedirs(out)
        cls.monomerize_obj = Monomerize(directory=os.path.join(cls.dir_input, 'examples_to_test5'), outdir=out)

    @classmethod
    def tearDownClass(cls):
        #[rmtree(a_path) for a_path in cls.test_path_list]
        pass

    def test_get_files(self):
        files = ['6epu.pdb', '6epv.pdb', '6epx.pdb', '6hi3.pdb']
        whole_files = sorted([os.path.join(self.align_obj.directory, x) for x in files])
        self.assertCountEqual(self.align_obj._get_files, whole_files)

    def test_get_mapfiles(self):
        self.assertCountEqual(self.align_obj._get_maplist, [])

    def test_monomerize_single(self):
        self.monomerize_obj.monomerize_single(file=os.path.join('tests', 'data_for_tests', 'examples_to_test5', 'Mpro-x0978.pdb'))
        testcases = ['Mpro-x2097_A.pdb', 'Mpro-x2097_A_event.cpp4']
        [self.assertTrue(os.path.exists(os.path.join('tests', 'data_for_tests', 'mono', x))) for x in testcases]

    def test_monomerize_all(self):
        self.monomerize_obj.monomerize_all()
        testcases = ['Mpro-x0978_A.pdb', 'Mpro-x0981_A.pdb',
                     'Mpro-x2097_A.pdb', 'Mpro-x2097_A_event.cpp4',
                     'Mpro-x2119_A.pdb', 'Mpro-x2119_A_event.cpp4']
        [self.assertTrue(os.path.exists(os.path.join('tests', 'data_for_tests', 'mono', x))) for x in testcases]

    def test_align_from_mono(self):
        a = Align(os.path.join('tests', 'data_for_tests', 'mono'), pdb_ref='', mono=True)
        a.align(out_dir=os.path.join('tests', 'data_for_tests', 'tmpa'))

        a_test_cases = ['Mpro-x0978_A_bound.pdb', 'Mpro-x0981_A_bound.pdb','Mpro-x2097_A_bound.pdb', 'Mpro-x2119_A_bound.pdb', 'Mpro-x2097_A_event.cpp4', 'Mpro-x2119_A_event.cpp4']
        a_test_exists = [os.path.exists(os.path.join('tests', 'data_for_tests', 'tmpa', x)) for x in a_test_cases]
        [self.assertTrue(x) for x in a_test_exists]

        a2 = Align(os.path.join('tests', 'data_for_tests', 'mono'), pdb_ref='Mpro-x2119', mono=True)
        a2.align(out_dir = os.path.join('tests', 'data_for_tests', 'tmpa2'))

        a2_test_cases = ['Mpro-x0978_A_bound.pdb', 'Mpro-x0981_A_bound.pdb','Mpro-x2097_A_bound.pdb', 'Mpro-x2119_A_bound.pdb', 'Mpro-x2097_A_event.cpp4', 'Mpro-x2119_A_event.cpp4']
        a2_test_exists = [os.path.exists(os.path.join('tests', 'data_for_tests', 'tmpa', x)) for x in a2_test_cases]
        [self.assertTrue(x) for x in a2_test_exists]

        # Also test if error state ruins things
        # error = Align(os.path.join('tests', 'data_for_tests', 'mono'), pdb_ref='Mpro-x2119', mono=False)

    def test_align_with_maps(self):
        self.align_obj_w_maps.align(out_dir=os.path.join('tests', 'data_for_tests', 'tmp_map'))
        map_test_cases = ['Mpro-x0978_bound.pdb', 'Mpro-x0981_bound.pdb','Mpro-x2097_bound.pdb', 'Mpro-x2119_bound.pdb', 'Mpro-x2097_event.cpp4', 'Mpro-x2119_event.cpp4']
        map_test_exists = [os.path.exists(os.path.join('tests', 'data_for_tests', 'tmp_map', x)) for x in map_test_cases]
        [self.assertTrue(x) for x in map_test_exists]


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
        #align_obj = Align(ATAD2_directory)
        #self.assertIsNotNone(align_obj._save_align)
        pass

    def test_aligns_are_correct_done(self):
        pass

if __name__ == '__main__':
    unittest.main()
