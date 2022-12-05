"""
A quick download functionality for the impatient.
"""

import warnings
from typing import Optional, Tuple, List, Dict, Any

import io
import os
import pandas as pd
import requests
import zipfile
from rdkit.Chem import PandasTools


class QuickDownloader:
    """
    This is simply a polished interface to the `api/download_structures` endpoint,
    which is the same as the download button on the Fragalysis website.
    Namely, it lacks the extended functionality of the `xcexporter` module.
    Instantiating it with the target name will download the zip file of these and can be interactive with via different
    methods.

    :cvar: fragalysis_api_url: The URL of the Fragalysis API.
    :cvar: api_data: The default options for the API endpoint. Can be also overridden in the constructor.
    :ivar: target_name: The name of the target as on the main Fragalysis page, e.g. 'Mpro', case sensitive
    :ivar: zf: The zip file object. see https://docs.python.org/3/library/zipfile.html

    The contents of the ``zipfile.ZipFile`` object stored in the attribute ``.zf`` can be written to disk
    with ``.write_all()``.

    The contents of a file within the zipfile can be accessed by subscripting the ``QuickDownloader`` instance,
    with a string that is contained in the filename, e.g. ``quick['metadata']``, will return the metadata.csv file
    content without having to waste time with filepaths.

    .. code-block:: python
        from fragalysis_api import QuickDownloader
        import pandas as pd
        from typing import List
        print(f'Default settings are: {QuickDownloader.api_data}')

        # Check if the target name is right
        target_names: List[str] = QuickDownloader.retrieve_target_names()
        target_name='Mpro'
        assert target_name in target_names, f'Target named "{target_name}" not found in the list of targets'

        # Download the data
        quick = QuickDownloader(target_name=target_name)
        quick.write_all(directory='downloads')
        hits: pd.DataFrame = quick.to_pandas(star_dummy=True)

        # Not all files have the reference pdb block, so if it does not the template is returned:
        reference_pdbblock: str = quick.reference_pdbblock

    The class method ``QuickDownloader.retrieve_target_data`` will download all the metadata for the targets.
    while ``QuickDownloader.retrieve_target_names`` will return their names.
    """
    fragalysis_api_url = 'fragalysis.diamond.ac.uk/api/download_structures/'
    api_data = {
        'proteins': '',
        'event_info': False,
        'sigmaa_info': False,
        'diff_info': False,
        'trans_matrix_info': False,
        'NAN': False,
        'mtz_info': False,
        'cif_info': False,
        'NAN2': False,
        'map_info': False,
        'single_sdf_file': True,
        'sdf_info': False,
        'pdb_info': False,
        'bound_info': True,
        'metadata_info': True,
        'smiles_info': True,
        'static_link': False,
        'file_url': ''}

    def __init__(self, target_name: str, **options):
        """
        Given a target download the zip file and store it in ``self.zf``.

        :param target_name: The name of the target as on the main Fragalysis page, e.g. 'Mpro', case sensitive
        :param options: Any of the options for the endpoint, e.g. ``event_info=True``, cf. ``cls.api_data``
        """
        self.target_name = target_name
        self.api_data = {options.get(k, v) for k, v in self.api_data.items()}
        url_response: requests.Response = requests.post(f'https://{self.fragalysis_api_url}',
                                                        json={'target_name': self.target_name, **self.api_data})
        url_response.raise_for_status()
        self.file_url: str = url_response.json()['file_url']
        response: requests.Response = requests.get(f"https://{self.fragalysis_api_url}?file_url={self.file_url}",
                                                   allow_redirects=True)
        response.raise_for_status()
        self.zf = zipfile.ZipFile(io.BytesIO(response.content), "r")

    def __getitem__(self, item: str) -> str:
        """
        Subscript via filename, e.g. ``quick['metadata']`` will return the metadata.csv file content.

        :param item: Part of the filename whose contents will be returned
        :return: The contents of the file whose name contains ``item``
        """
        for fileinfo in self.zf.infolist():
            if item in fileinfo.filename:
                return self.zf.read(fileinfo.filename).decode('utf8')
        else:
            raise KeyError(f'No file with {item} in the name found.')

    def __iter__(self) -> Tuple[str, str]:
        for fileinfo in self.zf.infolist():
            yield fileinfo.filename, self.zf.read(fileinfo.filename).decode('utf8')

    def __len__(self):
        """
        :return: The number of molecules-protein PDBs in the zip file.
        """
        return sum(['aligned/' in info.filename for info in self.zf.infolist()])

    def write_all(self, directory: Optional[str] = None):
        """
        Writes all the files within the zip file to disk in ``directory``.
        """
        if directory is None:
            directory = self.target_name
        if not os.path.exists(directory):
            os.makedirs(directory)
        for fileinfo in self.zf.infolist():
            if os.path.split(fileinfo.filename)[0] != '':
                os.makedirs(os.path.join(directory, os.path.split(fileinfo.filename)[0]), exist_ok=True)
            with open(os.path.join(directory, fileinfo.filename), 'w') as f:
                f.write(self.zf.read(fileinfo.filename).decode('utf8'))

    def to_pandas(self, star_dummy=True) -> pd.DataFrame:
        """
        Combine the metadata (``self['metadata]``) with sdf block (``self['combined.sdf']``),
        into a single pandas DataFrame.
        Peculiarly, Fragalysis stores dummy atoms as `Xe` instead of `*` in older SMILES, which is the standard.
        """
        # make a combined table
        # Fragalysis does not give attributes in the sdf entries. This is instead stored in metadata.csv.
        sdf_block = self['combined.sdf']
        df = PandasTools.LoadSDF(io.StringIO(sdf_block)).set_index('ID')
        try:
            metadata_block = self['metadata.csv'].replace('Xe', '*') if star_dummy else self['metadata.csv']
            df = pd.concat([df,
                            pd.read_csv(io.StringIO(metadata_block), index_col=0).set_index('crystal_name')
                            ], axis=1)
        except KeyError:
            warnings.warn('No metadata.csv found (legacy data). Returning only the sdf file.')
        return df

    @property
    def reference_pdbblock(self) -> str:
        """
        Not all files have the reference pdb block, so if it does not the template is returned.

        :return: The reference PDB for the target.
        """
        try:
            return self['reference']
        except KeyError:
            first_response: requests.Response = requests.get(f'https://{self.fragalysis_api_url}/api/targets/')
            first_response.raise_for_status()
            template_url = first_response.json()['results'][0]['template_protein']
            # /media/pdbs/ path:
            second_response: requests.Response = requests.get(f'https://{self.fragalysis_api_url}/{template_url}')
            second_response.raise_for_status()
            return second_response.text

    @classmethod
    def retrieve_target_data(cls) -> List[Dict[str, Any]]:
        """
        :return: A list of all the target metadata available on the Fragalysis API.
        """
        response: requests.Response = requests.get(f'https://{cls.fragalysis_api_url}/api/targets/')
        response.raise_for_status()
        return response.json()['results']

    @classmethod
    def retrieve_target_names(cls) -> List[str]:
        """
        :return: A list of all the target names available on the Fragalysis API.
        """
        return [target['title'] for target in cls.retrieve_target_data()]
