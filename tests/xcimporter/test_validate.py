import unittest
import os
from fragalysis_api import Validate, ValidatePDB


class ValidateTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        data_path = os.path.join('tests', 'data_for_tests')
        cls.test0_dir = os.path.join(data_path, 'examples_to_test0')
        cls.test1_dir = os.path.join(data_path, 'examples_to_test1')
        cls.test2_dir = os.path.join(data_path, 'examples_to_test2')


class CorrectInput(ValidateTest):

    @classmethod
    def setUpClass(cls):
        super(CorrectInput, cls).setUpClass()
        cls.correct_obj = Validate(cls.test0_dir)

    def test_if_pdbs_are_valid(self):
        """
        Test that there is no faults with this input (is a correct input)
        """
        self.assertTrue(self.correct_obj.is_pdbs_valid)
        self.assertTrue(not bool(self.correct_obj.validate_pdbs))

    def test_if_dir_exists(self):
        """
        Test that it knows that the dir exists
        """
        self.assertTrue(self.correct_obj.does_dir_exist)

    def test_if_there_is_a_pdb_in_dir(self):
        """
        Test it finds pdbs in the dir
        """
        self.assertTrue(self.correct_obj.is_there_a_pdb_in_dir)

    def test_finds_all_files_in_dir(self):
        """
        Test that if finds all the pdbs in the dir
        """
        self.assertCountEqual(self.correct_obj.get_files, [os.path.join(self.test0_dir, '6epv.pdb'),
                                                           os.path.join(self.test0_dir, '6hi3.pdb'),
                                                           os.path.join(self.test0_dir, '6epu.pdb'),
                                                           os.path.join(self.test0_dir, '6epx.pdb')])


class NoneExistingDir(ValidateTest):

    @classmethod
    def setUpClass(cls):
        super(NoneExistingDir, cls).setUpClass()
        cls.val_obj = Validate('not_a_dir')

    def test_pdbs_are_valid(self):
        """
        Pdbs should be valid even though there are none, as there are no invalid pdbs.
        """
        self.assertTrue(self.val_obj.is_pdbs_valid)
        self.assertTrue(not bool(self.val_obj.validate_pdbs))

    def test_if_dir_exists(self):
        """
        It should fail as the dir does not exist
        """
        self.assertFalse(self.val_obj.does_dir_exist)

    def test_if_there_is_a_pdb_in_dir(self):
        """
        It should fail as there are no pdbs in the dir that does not exist
        """
        self.assertFalse(self.val_obj.is_there_a_pdb_in_dir)


class ErrorInput(ValidateTest):

    @classmethod
    def setUpClass(cls):
        super(ErrorInput, cls).setUpClass()
        cls.val_obj = Validate(cls.test1_dir)

    def test_if_dir_exists(self):
        """
        Tests the dir exists
        """
        self.assertTrue(self.val_obj.does_dir_exist)

    def test_if_there_is_a_pdb_in_dir(self):
        """
        Tests the dir has pdbs
        """
        self.assertTrue(self.val_obj.is_there_a_pdb_in_dir)

    def test_finds_all_files_in_dir(self):
        """
        Test that if finds all the pdbs in the dir
        """
        self.assertCountEqual(self.val_obj.get_files,
                         [os.path.join(self.test1_dir, '5g1n (copy).pdb'),
                          os.path.join(self.test1_dir, '1.pdb'),
                          os.path.join(self.test1_dir, '5g1n.pdb'),
                          os.path.join(self.test1_dir, 'To_big_pdb.pdb'),
                          os.path.join(self.test1_dir, '5g1p (copy).pdb'),
                          os.path.join(self.test1_dir, '5g1o.pdb'),
                          os.path.join(self.test1_dir, 'Empty_pdb.pdb'),
                          os.path.join(self.test1_dir, '1234567890asdfghjklqwe.pdb'),
                          os.path.join(self.test1_dir, '5g1p.pdb'),
                          os.path.join(self.test1_dir, '5xzc.pdb')])

    def test_finds_invalid_pdbs(self):
        # Check it catches that there are invalid pdbs in the data set
        self.assertFalse(self.val_obj.is_pdbs_valid)
        # Ensure it caught them all
        self.assertCountEqual(self.val_obj.validate_pdbs,
                              ['5g1n (copy)', '1', 'To_big_pdb', '5g1p (copy)',
                               'Empty_pdb', '1234567890asdfghjklqwe'])

    def test_white_listing(self):
        """
        Test that it finds pdbs which contain characters which are not in the whitelist
        """
        self.assertFalse(
            ValidatePDB(os.path.join(self.test1_dir, '5g1n (copy).pdb')).does_pdb_name_contain_only_whitelist_char)
        self.assertFalse(
            ValidatePDB(os.path.join(self.test1_dir, '5g1p (copy).pdb')).does_pdb_name_contain_only_whitelist_char)

    def test_finds_pdbs_outside_allowed_size(self):
        """
        Test if it finds empty pdb files
        """
        self.assertFalse(
            ValidatePDB(os.path.join(self.test1_dir, 'Empty_pdb.pdb')).is_pdb_allowed_size)
        """
        Test if it finds pdbs larger than allowed size (5 mb)
        """
        self.assertFalse(
            ValidatePDB(os.path.join(self.test1_dir, 'To_big_pdb.pdb')).is_pdb_allowed_size)

    def test_finds_pdbs_with_allowed_name_sizes(self):
        """
        Test if it finds pdbs with a name with less than 4 characters
        """
        self.assertFalse(
            ValidatePDB(os.path.join(self.test1_dir, '1.pdb')).is_pdb_name_within_size_limit)
        """
        Test if it finds pdbs with a name larger than 20 characters
        """
        self.assertFalse(
            ValidatePDB(os.path.join(self.test1_dir, '1234567890asdfghjklqwe.pdb')).is_pdb_name_within_size_limit)


if __name__ == '__main__':
    unittest.main()