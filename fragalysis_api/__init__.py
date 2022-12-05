from .xcglobalscripts.set_config import ConfigSetup  # this requires `config.ini` in the working directory

# import into Fragalysis requires gemmi
try:
    from .xcimporter.validate import Validate, ValidatePDB
    from .xcimporter.conversion_pdb_mol import set_up, convert_small_AA_chains, copy_extra_files
    from .xcimporter.xc_utils import to_fragalysis_dir
    from .xcimporter.sites import Sites, contextualize_crystal_ligands
    from .xcimporter.xcimporter import xcimporter
    from .xcimporter.single_import import import_single_file
    from .xcimporter.align import Align
except ImportError:
    import warnings
    warnings.warn("The xcimporter functionality is unavailable without Gemmi", category=ImportWarning)

# extract from Fragalysis does not require gemmi
from .xcextracter.getdata import GetTargetsData, GetMoleculesData, GetPdbData, GetMolgroupData
from .xcextracter.frag_web_live import can_connect
from .xcextracter.xcextracter import xcextracter

from .xcanalyser.graphcreator import GraphRequest, xcgraphcreator
from .xcanalyser.xcanalyser import xcanalyser



