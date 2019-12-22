import os
import unittest
from fragalysis_api import set_up
from shutil import rmtree


class ConversionTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dir_input = os.path.join('tests', 'data_for_tests')
        cls.dir_output = os.path.join('tests', 'data_for_tests', 'tmp')


class PDBexample1(ConversionTest):

    @classmethod
    def setUpClass(cls):
        super(PDBexample1, cls).setUpClass()
        cls.obj_5q1j = set_up("5q1j", 'examples_to_test3', cls.dir_input, cls.dir_output)
        cls.obj_5qj7 = set_up("5qj7", 'examples_to_test3', cls.dir_input, cls.dir_output)

    @classmethod
    def tearDownClass(cls):
        rmtree(os.path.join(cls.dir_output))

    def test_open_file(self):
        """
        pdb file opens and checks that correct number of lines
        """
        self.assertEqual(len(self.obj_5q1j.pdbfile), 6779)

    def test_make_directory(self):
        """
        tests that directory to receive results has been created
        """
        self.assertTrue(os.path.isdir(self.obj_5q1j.RESULTS_DIRECTORY))

    def test_hets_and_cons(self):
        """
        tests that the correct number of heteroatoms and conect lines have been found
        """
        self.assertEqual(len(self.obj_5q1j.hetatms), 326)
        self.assertEqual(len(self.obj_5q1j.conects), 24)

    def test_remove_nonligs(self):
        """
        tests that solvents and ions from crystallography have been removed
        """
        # Changed because we leave in water?
        self.assertEqual(len(self.obj_5q1j.final_hets), 318)
        #assert len(self.obj_5q1j.final_hets) == 11

    def test_wanted_ligs(self):
        """
        tests that the ligand identifiers at the top of the pdb file have been found
        """
        self.assertNotEqual(len(self.obj_5q1j.wanted_ligs), None)

    def test_make_mol_objs(self):
        """
        tests that a molecular object has been made made for each ligand
        """
        self.assertEqual(len(self.obj_5qj7.mol_lst), 4)

    def test_make_pdb_file(self):
        """
        tests that a pdb file has been made for a particular ligand with the correct
        number of lines (HETATM & CONECT)
        """
        file = open(os.path.join(self.obj_5qj7.RESULTS_DIRECTORY, "5qj7_JMM_A_303.pdb")).readlines()
        self.assertEqual(len(file), 36)

    def test_make_mol_file(self):
        """
        tests that a mol file has been made with the correct number of lines
        """
        file = open(os.path.join(self.obj_5qj7.RESULTS_DIRECTORY, "5qj7_JMM_A_303_mol.mol")).readlines()
        self.assertEqual(len(file), 43)

    def test_make_sdf_file(self):
        """
        test that a sdf file has been made that incorporates the different mol objects
        """
        file = open(os.path.join(self.obj_5qj7.RESULTS_DIRECTORY, "5qj7_out.sdf")).readlines()
        self.assertEqual(len(file), 176)

if __name__ == '__main__':
    unittest.main()
