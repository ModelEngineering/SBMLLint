"""Access to SBMLLint configuration data.

Provides two functions.
  setConfiguration(<path>) - sets the value of the configuration dictionary.
  getConfiguration() - returns the current value of the configuration dictionary.
"""

import os
import yaml

import SBMLLint.common.constants as cn
from SBMLLint.common import msgs

_config_dict = None  # Configuration dictionary


def setConfiguration(path=cn.CFG_DEFAULT_PATH, fid=None):
  """
  :param str path: path to configuration file
  :param TextIOWrapper: overrides the path if present
  Updates _config_dict
  Notes:
    1. Changes string "True", "False" to booleans
    2. Inserts defaults for missing configuration keys
  """
  global _config_dict
  # Get the configuration file
  if fid is None:
    fid = open(path, "r")
  lines = fid.readlines()
  fid.close()
  lines = '\n'.join(lines)
  result = yaml.safe_load(lines)
  # Validate the section names
  for name in result.keys():
    if not name in cn.CFG_SECTIONS:
      msgs.error("Invalid section name: %s" % name)
  # Adjust values
  for k, v in result.items():
    if v == "True":
      result[k] = True
    if v == "False":
      result[k] = False
  for k, v in cn.CFG_DEFAULTS.items():
    if not k in result:
      result[k] = v 
  # Establish the configuration dictionary
  _config_dict = result

def getConfiguration():
  return _config_dict
#
setConfiguration(path=cn.CFG_DEFAULT_PATH)
