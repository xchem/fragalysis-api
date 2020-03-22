import unittest
from fragalysis_api.xcextracter import getdata


class GetPdbDataTest(unittest.TestCase):

    def test_GetPdbData_exists(self):
        search = getdata.GetPdbData()
        self.assertIsNotNone(search)

    def test_get_pdb_file(self):
        search = getdata.GetPdbData()
        self.assertIsNotNone(search.get_pdb_file('NUDT5A-x0114_1'))


if __name__ == '__main__':
    unittest.main()
