"""Reads yaml configuration file for SBMLLint."""
# TODO:
# 1. Provide standard interface to accessing configuration files
#    and placement of .smbllint_cfg
# 2. Incorporate implicits into structured_names

import yaml

with open(".sbmllint_cfg", "r") as fd:
  lines = fd.readlines()

lines = '\n'.join(lines)
result = yaml.safe_load(lines)
import pdb; pdb.set_trace()
