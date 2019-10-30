#!/usr/bin/env python

import json
import os

non_ligs = ['DMS', 'CL', 'MG', 'EDO', 'HOH', 'PEG', 'MLI', 'NI']
json.dump(non_ligs, open(os.path.join("../", "non_ligs"), "w"))
