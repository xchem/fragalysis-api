import unittest
import os
import align as AlignClass

relative_directory = 'data/ATAD2'
cwd = os.getcwd()
path = os.path.split(cwd)
directory = os.path.join(str(path[0]), relative_directory)


class Align_test(unittest.TestCase):
    align_obj = AlignClass.Align(directory)

    def test_get_files(self):
        files = ['6epu.pdb', '6epv.pdb', '6epx.pdb', '6hi3.pdb']
        whole_files = sorted([os.path.join(self.align_obj.directory, x) for x in files])
        self.assertEqual(self.align_obj.pdb_in_list, whole_files)





if __name__ == '__main__':
    unittest.main()
