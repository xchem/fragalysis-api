import glob
import os
import string

class Validate:
    def __init__(self, directory):
        self.directory = directory
        self.validate_pdbs = self.get_files

    @property
    def is_pdbs_valid(self):
        """
        Checks whether all pdbs are valid. If just one isn't it returns False

        :return: boolean answer
        """
        return not bool(self._fail_list)

    @property
    def validate_pdbs(self):
        """
        List of all pdbs which failed validation

        :return: list of all pdbs which failed validation
        """
        return self._fail_list

    @validate_pdbs.setter
    def validate_pdbs(self, pdb_file_list):

        fail_list = []
        for pdb in pdb_file_list:
            pdb_val = ValidatePDB(pdb)
            if not pdb_val.is_pdb_valid:
                fail_list.append(pdb_val.pdb_name)

        self._fail_list = fail_list

    @property
    def get_files(self):
        """Extracts a list of paths for all pdbs within the given directory.

        :return: list containing the path for each pdb in the given directory
        """
        return glob.glob(os.path.join(self.directory, "*.pdb"))

    @property
    def does_dir_exist(self):
        """
        Checks if the directory exists.

        :return: boolean answer
        """
        return bool(os.path.isdir(self.directory))

    @property
    def is_there_a_pdb_in_dir(self):
        """
        Checks if there is at least one pdb file in the directory.

        :return: boolean answer
        """
        return bool(self.get_files)


class ValidatePDB:

    def __init__(self, pdb_file):
        self.pdb_file = pdb_file
        self.pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]

    @property
    def is_pdb_valid(self):
        """
        If the pdb has failed just 1 check return False
        """
        return not (False in self._PDB_validations)

    @property
    def _PDB_validations(self):
        """
        Yields the results from the different types of validations.
        """
        yield self.is_pdb_allowed_size
        yield self.does_pdb_name_contain_only_whitelist_char
        yield self.is_pdb_name_within_size_limit
        # yield self.does_pdb_have_same_amount_of_protein_chains

    @property
    def is_pdb_allowed_size(self):
        """Checks whether a pdb file is larger than the allowed file size limit.

        :return: name of pdb if it failed to be smaller than the limit.
        """
        file_size_limit = 5.  # value in mb

        if os.path.getsize(self.pdb_file) <= 0:
            print(f'{self.pdb_name} file is to small.')
            return False

        if os.path.getsize(self.pdb_file) / 1000000 > file_size_limit:
            print(f'{self.pdb_name} file larger than limit of {file_size_limit}mb.')
            return False

        return True

    @property
    def does_pdb_name_contain_only_whitelist_char(self):
        """Checks whether a pdb files name contains characters which are not allowed.

        :return: name of pdb if it contains disallowed characters.
        """
        for character in self.pdb_name:
            if character not in string.ascii_letters+string.digits:
                print(f'{self.pdb_name} contains characters which are not allowed.\n'
                      f'Please only use ASCII or digits.')
                return False

        return True

    @property
    def is_pdb_name_within_size_limit(self):
        """Checks whether the name of a pdb file is within the allowed range.

        :return: name of pdb if it was outside the allowed range.
        """
        char_min = 4
        char_max = 20

        if (len(self.pdb_name) > char_max) or (len(self.pdb_name) < char_min):
            print(f'{self.pdb_name} has a name outside the allowed length range.\n'
                  f'Please keep names between {char_min}-{char_max} characters.')
            return False

        return True
