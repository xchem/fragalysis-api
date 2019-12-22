import os
import unittest
from fragalysis_api import to_fragalysis_dir
from glob import glob
from shutil import rmtree


class xcUtilsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dir_input = os.path.join('tests', 'data_for_tests')


class ToFragalysisDir(xcUtilsTest):

    @classmethod
    def setUpClass(cls):
        super(ToFragalysisDir, cls).setUpClass()
        to_fragalysis_dir('a_test', os.path.join(cls.dir_input, 'examples_to_test4'))

    @classmethod
    def tearDownClass(cls):
        a_test_path_list = glob(os.path.join(cls.dir_input, '*a_test*'))
        [rmtree(a_path) for a_path in a_test_path_list]


    def test_hey(self):
        pass
