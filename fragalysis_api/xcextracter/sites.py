from fragalysis_api import GetMoleculesData, GetTargetsData, ConfigSetup
import pandas as pd
import urllib
import json


class GetMolgroupData:
    '''
    Class to contain and get data about molecules in  a fragalysis site
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
        self.molgroup_url = settings.get('molgroup_from_description', 'search')
        self.query = settings.get('molgroup_from_description', 'query')

        # get full url
        self.search_url = str(self.frag_url + self.molgroup_url + self.query)
        self.target_name = None
        self.target_name_url = None
        self.target_json = None
        self.target_id_list = None
        self.molgroup_name_url = None
        self.mol_ids = None
        self.site_mols = None

    def set_molgroup_url(self, target, description):
        '''
        Setting target name url

        :param target: Target name
        :param description: Site name
        '''
        self.target_name = target
        search_description = description.replace(' ', '+')
        self.molgroup_name_url = str(self.search_url + search_description)

        target_search = GetTargetsData()
        target_search.set_target_name_url(target=target)

        target_search.get_target_json()
        ids = target_search.get_target_id_list()
        print(ids)
        if len(ids) != 1:
            raise Exception('More than one target found with that name')
        else:
            target_id = ids[0]
            self.molgroup_name_url = self.molgroup_name_url + '&target_id=' + str(target_id)

    def get_molgroup_json(self):
        if not self.molgroup_name_url:
            raise Exception('Please initiate target_name url with set_molgroup_url(<description>)!')

        # get response from url and decode -> json
        with urllib.request.urlopen(self.molgroup_name_url) as f:
            response = json.loads(f.read().decode('utf-8'))

        # set json as decoded response for processing
        self.target_json = response

    def extract_molids(self):
        if len(self.target_json['results']) != 1:
            raise Exception('Either the query has not been run, or too many results were returned')
        else:
            mol_ids = self.target_json['results'][0]['mol_id']
            self.mol_ids = mol_ids

    def get_site_molinfo(self):
        search = GetMoleculesData()
        search.set_target_id(target=self.target_name)
        search.set_molecule_url()
        search.set_mol_data()
        search.set_complete_mol_data()
        data = search.get_complete_mol_data
        data2 = [x for x in data if x['id'] in self.mol_ids]
        self.site_mols = data2

    def get_site_mol_table(self, target, site_name):
        self.set_molgroup_url(target=target, description=site_name)
        self.get_molgroup_json()
        self.extract_molids()
        self.get_site_molinfo()
        return pd.DataFrame(self.site_mols)
