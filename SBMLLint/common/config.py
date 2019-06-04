
CONFIG_PATH = "SBMLLint/.sbmllint_cfg"
config_dict = {}  # Default configuration dictionary.

"""Reads yaml configuration file for SBMLLint."""
# TODO:
# 1. Provide standard interface to accessing configuration files
#    and placement of .smbllint_cfg
# 2. Incorporate implicits into structured_names

import yaml

import SBMLLint.common.constants as cn

def getConfiguration(path=None):
  """
  :param str path: path to configuration file
  :return dict: dictionary of configuration values
  Notes:
    1. Changes string "True", "False" to booleans
    2. Inserts defaults for missing configuration keys
  """
  # TODO: Consider that .sbmllint_cfg may be in home directory
  if path is None:
    path = CONFIG_PATH
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
  return result
#
config_dict = getConfiguration()
