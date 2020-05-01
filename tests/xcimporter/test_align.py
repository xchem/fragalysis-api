import unittest
import os
from fragalysis_api import Align
from glob import glob
from shutil import rmtree


class AlignTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dir_input = os.path.join('tests', 'data_for_tests')


class EasyAlign(AlignTest):

    @classmethod
    def setUpClass(cls):
        super(EasyAlign, cls).setUpClass()
        cls.align_obj = Align(os.path.join(cls.dir_input, 'examples_to_test0'), pdb_ref='')
        cls.align_obj_w_ref = Align(os.path.join(cls.dir_input, 'examples_to_test0'), pdb_ref='6hi3')
        cls.align_obj_w_wrong_ref = Align(os.path.join(cls.dir_input, 'examples_to_test0'), pdb_ref='wrong_pdb')

    @classmethod
    def tearDownClass(cls):
        #[rmtree(a_path) for a_path in cls.test_path_list]
        pass

    def test_get_files(self):
        files = ['6epu.pdb', '6epv.pdb', '6epx.pdb', '6hi3.pdb']
        whole_files = sorted([os.path.join(self.align_obj.directory, x) for x in files])
        self.assertCountEqual(self.align_obj._get_files, whole_files)

    def test_load_obj_successfully(self):
        self.assertIsNotNone(self.align_obj._load_objs)

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
