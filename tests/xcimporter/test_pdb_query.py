import unittest
from fragalysis_api.xcimporter.pdbimporter import ImportPdb
from fragalysis_api.xcimporter.pdbquery import Query
import os


class ValidateTest(unittest.TestCase):

    def test_query(self):
        print('testing query')
        pdb_code = '6epu'
        chain_id = 'A'
        query_obj = Query(pdb_code, chain_id)
        query_obj.get_matching_proteins()
        query_obj.get_ligands()
        assert(len(query_obj.match_ligs) == 41)

if __name__ == '__main__':
    unittest.main()