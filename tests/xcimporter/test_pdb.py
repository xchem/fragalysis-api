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
        self.assertTrue(obj.user_id is not None)
        self.assertTrue(obj.pdb_code is not None)
        self.assertTrue(obj.data_dir is not None)
        #obj.pdb_importer()
        self.assertFalse(obj.pdb_exists)

if __name__ == '__main__':
    unittest.main()

