import unittest
import os
from fragalysis_api import Validate, ValidatePDB

data_path = os.path.join('data', 'xcimporter', 'input')
good_dir = os.path.join(data_path, 'ATAD2')
test_dir = os.path.join(data_path, 'examples_to_test1')
test2_dir = os.path.join(data_path, 'examples_to_test2')
PDB_dir = os.path.join(data_path, 'PDB')


class ValidateTest(unittest.TestCase):

    def test_good_input(self):
        """
        This test shows how it should work when the input is correct
        """
        val_obj = Validate(good_dir)
        # Test that there is no faults with this input
        self.assertTrue(val_obj.is_pdbs_valid)
        self.assertTrue(not bool(val_obj.validate_pdbs))
        # If dir exists
        self.assertTrue(val_obj.does_dir_exist)
        # If there is a pdb in dir
        self.assertTrue(val_obj.is_there_a_pdb_in_dir)
        # Does it get all the files?
        self.assertCountEqual(val_obj.get_files, [os.path.join(good_dir, '6epv.pdb'), os.path.join(good_dir, '6hi3.pdb'),
                                             os.path.join(good_dir, '6epu.pdb'),  os.path.join(good_dir, '6epx.pdb')])


    def test_miscellanea(self):
        """
        When the input is a cif file + no pdbs in the directory
        """
        misc_obj = Validate('not_a_dir')
        # Check pdbs are valid (there are no pdbs so it should pass)
        self.assertTrue(misc_obj.is_pdbs_valid)
        self.assertTrue(not bool(misc_obj.validate_pdbs))
        # If dir exists
        self.assertFalse(misc_obj.does_dir_exist)
        # That there are no pdbs in dir
        self.assertFalse(misc_obj.is_there_a_pdb_in_dir)

    def test_mistake_catching(self):
        """
        Most of these cases should raise errors
        """
        hard_obj = Validate(test_dir)
        # If dir exists
        self.assertTrue(hard_obj.does_dir_exist)
        # If there are a pdb in dir
        self.assertTrue(hard_obj.is_there_a_pdb_in_dir)
        # Check all files are retrieved
        self.assertCountEqual(hard_obj.get_files,
                         [os.path.join(test_dir, '5g1n (copy).pdb'), os.path.join(test_dir, '1.pdb'),
                          os.path.join(test_dir, '5g1n.pdb'), os.path.join(test_dir, 'To_big_pdb.pdb'),
                          os.path.join(test_dir, '5g1p (copy).pdb'), os.path.join(test_dir, '5g1o.pdb'),
                          os.path.join(test_dir, 'Empty_pdb.pdb'), os.path.join(test_dir, '1234567890asdfghjklqwe.pdb'),
                          os.path.join(test_dir, '5g1p.pdb'), os.path.join(test_dir, '5xzc.pdb')])
        # Check it catches that there are invalid pdbs in the data set
        self.assertFalse(hard_obj.is_pdbs_valid)
        # Ensure it caught them all
        self.assertCountEqual(hard_obj.validate_pdbs, ['5g1n (copy)', '1', 'To_big_pdb', '5g1p (copy)',
                                                  'Empty_pdb', '1234567890asdfghjklqwe'])
        # Check if whitelisting works
        self.assertFalse(
            ValidatePDB(os.path.join(test_dir, '5g1n (copy).pdb')).does_pdb_name_contain_only_whitelist_char)
        self.assertFalse(
            ValidatePDB(os.path.join(test_dir, '5g1p (copy).pdb')).does_pdb_name_contain_only_whitelist_char)
        # Check if it finds empty pdb files
        self.assertFalse(
            ValidatePDB(os.path.join(test_dir, 'Empty_pdb.pdb')).is_pdb_allowed_size)
        # Check if it finds to large files
        self.assertFalse(
            ValidatePDB(os.path.join(test_dir, 'To_big_pdb.pdb')).is_pdb_allowed_size)
        # Check if it finds pdbs with too small names
        self.assertFalse(
            ValidatePDB(os.path.join(test_dir, '1.pdb')).is_pdb_name_within_size_limit)
        # Check if it finds pdbs with too large names
        self.assertFalse(
            ValidatePDB(os.path.join(test_dir, '1234567890asdfghjklqwe.pdb')).is_pdb_name_within_size_limit)

    def test_semi_hard(self):
        pass

if __name__ == '__main__':
    unittest.main()