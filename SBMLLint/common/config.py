"""Access to SBMLLint configuration data.

Provides two functions.
  setConfiguration(<path>) - sets the value of the configuration dictionary.
  getConfiguration() - returns the current value of the configuration dictionary.
"""

import os
import yaml

import SBMLLint.common.constants as cn

_config_dict = None  # Configuration dictionary


def setConfiguration(path=cn.CFG_DEFAULT_PATH):
  """
  :param str path: path to configuration file
  Updates _config_dict
  Notes:
    1. Changes string "True", "False" to booleans
    2. Inserts defaults for missing configuration keys
  """
  # TODO: Consider that .sbmllint_cfg may be in home directory
  global _config_dict
  with open(path, "r") as fd:
    lines = fd.readlines()
  lines = '\n'.join(lines)
  result = yaml.safe_load(lines)
  for k, v in result.items():
    if v == "True":
      result[k] = True
    if v == "False":
      result[k] = False
  for k, v in cn.CFG_DEFAULTS.items():
    if not k in result:
      result[k] = v 
  _config_dict = result

def getConfiguration():
  return _config_dict
#
setConfiguration(path=cn.CFG_DEFAULT_PATH)
