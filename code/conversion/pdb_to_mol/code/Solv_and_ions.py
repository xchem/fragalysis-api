#!/usr/bin/env python

import json
import os

non_ligs = ['DMS', 'CL', 'MG', 'EDO', 'HOH', 'PEG']
json.dump(non_ligs, open(os.path.join("./conversion", "non_ligs"), "w"))
