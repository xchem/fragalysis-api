import unittest
from fragalysis_api.xcimporter.pdbimporter import ImportPdb
from fragalysis_api.xcimporter.pdbquery import Query
import os


class ValidateTest(unittest.TestCase):

    def test_download(self):
        data_dir = os.path.join('data', 'xcimporter', 'input')
        user_id = 'examples_to_test3'
        pdb_code = '6epu'
        if os.path.exists(os.path.join(data_dir, user_id, pdb_code+'.pdb')):
            os.remove(os.path.join(data_dir, user_id, pdb_code+'.pdb'))
        obj = ImportPdb(data_dir, user_id, pdb_code)
        self.assertTrue(obj.pdb_importer())

    def test_query(self):
        in_dir = os.path.join('data', 'xcimporter', 'input')
        out_dir = os.path.join('data', 'xcimporter', 'output')
        pdb_code = '6epu'
        chain_id = 'A'
        user_id = 'examples_to_test3'
        query_obj = Query(pdb_code, chain_id)
        query_obj.get_ligands()
        assert(len(query_obj.match_ligs) == 41)


if __name__ == '__main__':
    unittest.main()

