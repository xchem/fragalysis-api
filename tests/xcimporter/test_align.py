import unittest
import os
from fragalysis_api import Align

# Load test data

cwd = os.getcwd()
path = os.path.split(cwd)

data_path = os.path.join(str(path[0]), '..', 'data', 'xcimporter', 'input')
ATAD2_directory = os.path.join(data_path, 'ATAD2')
cif_directory = os.path.join(data_path, 'input/CIF')
Hard_directory = os.path.join(data_path, 'Hard_example')
Semi_hard_directory = os.path.join(data_path, 'Semi_hard_examples')
PDB_directory = os.path.join(data_path, 'PDB')

directory = os.path.join(str(path[0]), ATAD2_directory)


class Align_test(unittest.TestCase):


    # def test_get_files(self):
    #     files = ['6epu.pdb', '6epv.pdb', '6epx.pdb', '6hi3.pdb']
    #     whole_files = sorted([os.path.join(self.align_obj.directory, x) for x in files])
    #     self.assertEqual(self.align_obj.pdb_in_list, whole_files)

    def test_load_success(self):
        align_obj = Align(ATAD2_directory)
        self.assertIsNotNone(align_obj._load_objs)

    def test_get_ref_success(self):
        align_obj = Align(ATAD2_directory)
        self.assertTrue(align_obj._get_ref, '6epv.pdb')

    def test_align_success(self):
        align_obj = Align(ATAD2_directory)
        self.assertIsNotNone(align_obj._save_align)



if __name__ == '__main__':
    unittest.main()
