
CONFIG_PATH = "SBMLLint/.sbmllint_cfg"

"""Reads yaml configuration file for SBMLLint."""
# TODO:
# 1. Provide standard interface to accessing configuration files
#    and placement of .smbllint_cfg
# 2. Incorporate implicits into structured_names

import yaml

def getConfiguration(path=None):
  """
  :param str path: path to configuration file
  :return dict: dictionary of configuration values
  """
  # TODO: Consider that .sbmllint_cfg may be in home directory
  if path is None:
    path = CONFIG_PATH
  with open(path, "r") as fd:
    lines = fd.readlines()
  lines = '\n'.join(lines)
  result = yaml.safe_load(lines)
  return result
