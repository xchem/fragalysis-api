import glob
import os
import string


class Validate:
    def __init__(self, directory):
        self.directory = directory
        self.validate_pdbs = self.get_files

    def is_pdbs_val(self):
        return if 

    @property
    def validate_pdbs(self):

        return self._fail_list

    @validate_pdbs.setter
    def validate_pdbs(self, pdb_file_list):

        fail_list = []
        for pdb in pdb_file_list:
            self.is_directory_empty(pdb)
            fail_list += [i for i in ValidatePDB(pdb).PDB_validations() if i is not None]

        self._fail_list = fail_list

    @property
    def get_files(self):
        """Extracts a list of paths for all pdbs within the given directory.

        :return: list containing the path for each pdb in the given directory
        """
        return glob.glob(os.path.join(self.directory, "*.pdb"))

    def is_directory_empty(self, pdb_file):
        pass


class ValidatePDB:

    def __init__(self, pdb_file):
        self.pdb_file = pdb_file
        self.pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]


    def PDB_validations(self):
        yield self.is_pdb_to_large()
        yield self.does_pdb_name_contain_none_whitelist_char()
        yield self.is_pdb_name_within_range()
        yield self.does_pdb_have_same_amount_of_protein_chains()


    def is_pdb_to_large(self):
        """Checks whether a pdb file is larger than the allowed file size limit.

        :return: name of pdb if it failed to be smaller than the limit.
        """
        file_size_limit = 5.  # value in mb

        if os.path.getsize(self.pdb_file) <= 0:
            print(f'{self.pdb_name} file is to small.')
            return self.pdb_name

        if os.path.getsize(self.pdb_file) / 1000000 > file_size_limit:
            print(f'{self.pdb_name} file larger than limit of {file_size_limit}mb.')
            return self.pdb_name

    def does_pdb_name_contain_none_whitelist_char(self):
        """Checks whether a pdb files name contains characters which are not allowed.

        :return: name of pdb if it contains disallowed characters.
        """
        for character in self.pdb_name:
            if character not in string.ascii_letters+string.digits:
                print(f'{self.pdb_name} contains characters which are not allowed.\n'
                      f'Please only use ASCII or digits.')
                return self.pdb_name

    def is_pdb_name_within_range(self):
        """Checks whether the name of a pdb file is within the allowed range.

        :return: name of pdb if it was outside the allowed range.
        """
        char_min = 4
        char_max = 20

        if (len(self.pdb_name) > char_max) or (len(self.pdb_name) < char_min):
            print(f'{self.pdb_name} has a name outside the allowed length range.\n'
                  f'Please keep names between {char_min}-{char_max} characters.')
            return self.pdb_name

    def does_pdb_have_same_amount_of_protein_chains(self):
        pass



#if __name__ == '__main__':

 #   a_dir = '../anna/Hard_example/'
  #  validation = Validate(a_dir)
   # print(validation.get_files)
    #print(validation.validate_pdbs)
