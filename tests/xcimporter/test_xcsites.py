import unittest
import os
from fragalysis_api import xcimporter, Sites, contextualize_crystal_ligands


class XcSitesTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.in_dir = os.path.join(
            'tests', 'data_for_tests', 'examples_to_test5')
        cls.in_file = os.path.join(
            'tests', 'data_for_tests', 'examples_to_test5', 'Mpro-x0978.pdb')
        cls.in_file_no_ext = 'Mpro-x0978'
        cls.out_dir = os.path.join('tests', 'data_for_tests')
        cls.target = 'TestTarget'
        cls.validate = False
        cls.rrf = True
        cls.biomol = os.path.join(
            'tests', 'data_for_tests', 'examples_to_test5', 'biomol.txt')
        cls.metadata = True
        cls.covalent = True
        cls.reference = 'Mpro-x0981_A'
        cls.mll = -1

    @classmethod
    def tearDownClass(cls):
        #rmtree(os.path.join(cls.out_dir, cls.target))
        pass

    def test_sites(self):

        xcimporter(
            in_dir=self.in_dir,
            out_dir=self.out_dir,
            target=self.target,
            validate=self.validate,
            reduce_reference_frame=self.rrf,
            metadata=self.metadata,
            biomol=self.biomol,
            covalent=self.covalent,
            pdb_ref=self.reference,
            max_lig_len=self.mll
        )

        folder = os.path.join(self.out_dir, self.target)
        site_obj = Sites.from_folder(folder, recalculate=False)
        # site_obj.cluster_missing_mols(
        #        folder=folder, com_tolerance=5.00, other_tolerance=1.00)
        site_obj.to_json()
        contextualize_crystal_ligands(folder=folder)
        site_obj.apply_to_metadata()

        self.assertTrue(expr=os.path.exists(os.path.join(
            self.out_dir, self.target, 'aligned', 'Mpro-x0978_0A', 'Mpro-x0978_0A_sites.json')))


if __name__ == '__main__':
    unittest.main()
