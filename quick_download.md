## Quick download

This is not part of the regular package, which allows more precise control, but simply a wrapper to
the `api/download_structures` endpoint (which is used by the download modal) to download all the data
for a given target.

This endpoint differs from the links in the `api/targets` endpoint response. Namely,
for say NUDT7A, the JSON response to https://fragalysis.diamond.ac.uk/api/targets/?title=NUDT7A
will feature `zip_archive` that will point to `https://fragalysis.diamond.ac.uk/media/targets/👾👾👾.zip`
and a `template_protein`, which is (possibly) the reference structure for molecular replacement and 
not necessary for alignment. NUDT7 is an older target and actually lacks a reference structure in the download.

```python
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
```
