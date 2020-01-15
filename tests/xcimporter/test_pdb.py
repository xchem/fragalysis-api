import unittest
from fragalysis_api.xcimporter.pdbimporter import ImportPdb
import os


class ValidateTest(unittest.TestCase):

    def test_download(self):
        print('testing download')
        data_dir = os.path.join('data', 'xcimporter', 'input')
        user_id = 'examples_to_test3'
        pdb_code = '6epu'
        if os.path.exists(os.path.join(data_dir, user_id, pdb_code+'.pdb')):
            os.remove(os.path.join(data_dir, user_id, pdb_code+'.pdb'))
        obj = ImportPdb(data_dir, user_id, pdb_code)
        self.assertFalse(obj.user_id is None)
        self.assertFalse(obj.pdb_code is None)
        self.assertFalse(obj.data_dir is None)
       # self.assertTrue(obj.pdb_importer())

if __name__ == '__main__':
    unittest.main()

