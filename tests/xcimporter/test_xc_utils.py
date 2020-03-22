import os
import unittest
from fragalysis_api import to_fragalysis_dir
from glob import glob
from shutil import rmtree


class xcUtilsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dir_input = os.path.join('tests', 'data_for_tests')


class ToFragalysisDir(xcUtilsTest):

    @classmethod
    def setUpClass(cls):
        super(ToFragalysisDir, cls).setUpClass()
        to_fragalysis_dir('a_test', os.path.join(cls.dir_input, 'examples_to_test4'))
        cls.test_path_list = glob(os.path.join(cls.dir_input, '*a_test*'))

    @classmethod
    def tearDownClass(cls):
        [rmtree(a_path) for a_path in cls.test_path_list]

    def test_correct_dirs_created(self):
        """
        Tests that the correct sub-dirs gets created
        """
        self.assertCountEqual(self.test_path_list, [os.path.join(self.dir_input, 'a_test_6epx'),
                                                    os.path.join(self.dir_input, 'a_test_6epv'),
                                                    os.path.join(self.dir_input, 'a_test_6hi3'),
                                                    os.path.join(self.dir_input, 'a_test_6epu')])

    def test_only_files_from_the_same_pdb_in_subdirs(self):
        """
        Tests each dir only has files from the same pdb
        """
        pass

    def test_all_files_are_allocated_to_a_new_dir(self):
        pass

    def test_apo_files_gets_created(self):

        pass

    def test_bound_files_gets_created(self):
        pass

    def test_mol_and_sdf_files_gets_created(self):
        pass
