import unittest
import os
from fragalysis_api import Validate, ValidatePDB

cwd = os.getcwd()
path = os.path.split(cwd)

data_path = os.path.join('data', 'xcimporter', 'input') #str(path[0]), '..', 
good_directory = os.path.join(data_path, 'ATAD2')
cif_directory = os.path.join(data_path, 'CIF')
Hard_directory = os.path.join(data_path, 'Hard_example')
Semi_hard_directory = os.path.join(data_path, 'Semi_hard_examples')
PDB_directory = os.path.join(data_path, 'PDB')


class Validate_test(unittest.TestCase):

    def test_good_input(self):
        """
        This test shows how it should work when the input is correct
        """
        print
        val_obj = Validate(good_directory)
        # Test that there is no faults with this input
        self.assertTrue(val_obj.is_pdbs_valid)
        self.assertTrue(not bool(val_obj.validate_pdbs))
        # Test that the dir contains files (it doesnt return a mistake)
        self.assertTrue(val_obj.does_dir_exist)
        self.assertTrue(val_obj.is_there_a_pdb_in_dir)

    def test_cif(self):
        """
        When the input is a cif file + no pdbs in the directory
        """
        cif_obj = Validate(cif_directory)
        # Test that there is no faults with this input
        self.assertTrue(cif_obj.is_pdbs_valid)
        self.assertTrue(not bool(cif_obj.validate_pdbs))
        # Test that the dir contains files (it doesnt return a mistake)
        self.assertTrue(cif_obj.does_dir_exist)
        self.assertFalse(cif_obj.is_there_a_pdb_in_dir)


if __name__ == '__main__':
    unittest.main()