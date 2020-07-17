import urllib
import json
from fragalysis_api import ConfigSetup


class GetTargetsData:
    '''
    Class to contain and get data on the protein target of interest
    '''
    def __init__(self):
        '''
        :param self.frag_url: URL of the fragalysis website
        :param self.target_url: URL of extention for the target of fragalysis website
        :param self.query: URL to tell the restfull API how to do the query
        :param self.search_url: URL for the query 
        :param self.target_name_url: Lists of all the proteins associated with the target
        :param self.target_json: Json of the query
        :param self.target_id_list: list of the ID numbers/number for the targets/target
        '''
        settings = ConfigSetup()

        self.frag_url = settings.get('fragalysis', 'url')
        self.target_url = settings.get('targets', 'search')
        self.query = settings.get('targets', 'query')

        # get full url
        self.search_url = str(self.frag_url + self.target_url + self.query)

        self.target_name_url = None
        self.target_json = None
        self.target_id_list = None

    def set_target_name_url(self, target):
        '''
        Setting target name url

        :param target: Target name
        '''
        self.target_name_url = str(self.search_url + target)

    def get_target_json(self):
        '''
        Gets the json output of the lists of all the proteins associated with the target
        from fragalysis

        :returns response: Json response from the restful API
        '''
        if not self.target_name_url:
            raise Exception('Please initiate target_name url with set_target_name_url(<target>)!')

        # get response from url and decode -> json
        with urllib.request.urlopen(self.target_name_url) as f:
            response = json.loads(f.read().decode('utf-8'))

        # set json as decoded response for processing
        self.target_json = response

        return response

    def get_target_id_list(self):
        '''
        Gets the target_id_ from self.target_json associated with the queried target in the 
        fragalysis database

        :returns id_list: The ID number associated with the target.
        '''
        if not self.target_json:
            raise Exception('Please get data with get_target_json!')

        id_list = [result['id'] for result in self.target_json['results']]
        self.target_id_list = id_list
        return id_list


class GetPdbData:
    '''
    Getting protine data from the fragalysis website in the pdb file format
    '''
    def __init__(self):
        '''
        :param self.frag_url: URL of the fragalysis website
        :param self.pdb_url: URL to search the fragalysis data base for the protein
        :param self.query: URL to tell the restfull API how to do the query
        '''
        settings = ConfigSetup()

        self.frag_url = settings.get('fragalysis', 'url')
        self.pdb_url = settings.get('pdb', 'search')
        self.bound_pdb_url = settings.get('bound_pdb', 'search')
        self.bound_query = settings.get('bound_pdb', 'query')
        self.query = settings.get('pdb', 'query')

    def get_apo_pdb_file(self, code):
        '''
        Function to search fragalysis for PDB file.
        :param code: The code associated with the protein on the fragalysis website.
        :return response: Response from the API = pdb file
        '''
        url = str(self.frag_url + self.pdb_url + self.query + code)
        # get response from url and decode -> json
        with urllib.request.urlopen(url) as f:
            response = json.loads(f.read().decode('utf-8'))

        if len(response['results']) > 1:
            raise Exception('more than one pdb found... contact admin')
        return response['results'][0]['pdb_data']

    def get_bound_pdb_file(self, code):
        '''
        Function to search fragalysis for PDB file.
        :param code: The code associated with the protein on the fragalysis website.
        :return response: Response from the API = pdb file
        '''
        url = str(self.frag_url + self.bound_pdb_url + self.bound_query + code)
        # get response from url and decode -> json
        with urllib.request.urlopen(url) as f:
            response = json.loads(f.read().decode('utf-8'))

        if len(response['results']) > 1:
            raise Exception('more than one pdb found... contact admin')
        return response['results'][0]['bound_pdb_data']


class GetMoleculesData:
    def __init__(self):
        '''
        :param self.frag_url: URL of the fragalysis website
        :param self.moleclues_url: URL of extention for the moleclue on the fragalysis website
        :param self.query: URL to tell the restfull API how to do the query
        :param self.search_url: URL for the query 
        :param self.get_molecule_url: Molecule URL for the fragalysis website
        :param self.get_target_id: Target's fragalysis ID 
        :param self.get_mol_data: Molecule data in json format
        :param self.get_complete_mol_data: Complete molecule data in json format
        '''
        settings = ConfigSetup()


        self.frag_url = settings.get('fragalysis', 'url')
        self.molecules_url = settings.get('molecules', 'search')
        self.query = settings.get('molecules', 'query')

        self.search_url = str(self.frag_url + self.molecules_url + self.query)

        self.get_molecule_url = ''
        self.get_target_id = ''
        self.get_mol_data = ''
        self.get_complete_mol_data = ''

    def set_target_id(self, target):
        """
        Gets the targets fragalysis ID (a number) based on the target
        :param target: name of target (e.g. ATAD)
        :return: initializes the target ID in self.taget_ids
        """
        search = GetTargetsData()
        search.set_target_name_url(target)
        search.get_target_json()
        id_list = search.get_target_id_list()
        self.get_target_id = id_list[0]  # Change id_list to be a single id (will never be more than 1)

    def set_molecule_url(self):
        """
        Sets the molecule url based on the targets fragalysis ID
        :return:
        """
        url = str(self.search_url + str(self.get_target_id))

        self.get_molecule_url = url

    def set_mol_data(self):
        """
        Sets the molecule data in json format
        :return:
        """

        # get response from url and decode -> json
        print(self.get_molecule_url)
        with urllib.request.urlopen(self.get_molecule_url) as f:
            self.get_mol_data = json.loads(f.read().decode('utf-8'))

    def set_complete_mol_data(self):
        '''
        Extracting the results from the url query
        i.e self.get_complete_mol_data contains only
        the results of the query on te molecule.
        :return:
        '''

        if not self.get_target_id:
            raise Exception('Please get the target ids with get_target_ids!')

        json_list = []

        json_list.extend(self.get_mol_data['results'])

        self.get_complete_mol_data = json_list
        