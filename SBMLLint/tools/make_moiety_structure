#!/usr/bin/env python
"""
Creates the moiety_structure section of the SBMLLint
configuration file. The user declares all moieties.
The moiety structure for a molecule is the number of
occurrences of each moiety string in the molecule.
Substrings are handled by looking first for the longest
moiety names.
"""

import argparse
import yaml
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common.moiety import MoietyStoichiometry
from SBMLLint.common import msgs


DEFAULT_CONFIG_FILE = "config.yml"
MOIETY_COUNT_SEPARATOR = ", "


def getMoieties(moiety_fid):
  """
  Creates the YAML configuration file.
  :param IOWrappter moiety_fid: 
  :return list-str: list of moiety names
  :raises ValueError: if one moiety is a substring of another
  """
  return yaml.safe_load(moiety_fid)

# TODO: Extend to do recursive checks of moiety names for
# their substructures
def findMoietyStoichiometries(molecule_name, moiety_names,
    is_check_error=True):
  """
  Finds the moieties and their stoichiometries for each molecule.
  Handles substrings by looking for largest strings first.
  :param str molecule:
  :param list-str moiety_names:
  :return list-MoietyStoichiometry:
  """
  result = []
  sorted_moiety_names = sorted(
      moiety_names, key=lambda n: len(str(n)), reverse=True)
  current_name = molecule_name
  for moiety_name in sorted_moiety_names:
    moiety_name = str(moiety_name)
    count = current_name.count(moiety_name)
    if count > 0:
      current_name = current_name.replace(moiety_name, "")
      result.append(MoietyStoichiometry(moiety_name, count,
          MOIETY_COUNT_SEPARATOR))
  current_name = current_name.replace("_", "")
  if is_check_error and len(current_name) > 0:
    msg = "Moieties specified do not completely express %s" % molecule_name
    msgs.error(msg)
  return result

# TODO: Look for moieties from largest string to smallest,
# eliminating those found.
def main(xml_fid, moiety_fid, config_fid):
  """
  Creates the YAML configuration file.
  :param IOWrapper xml_fid: xml file:
  :param IOWrapper moiety_fid: 
  :param IOWrapper config_fid: file to write
  """
  # Acquire and validate moieties
  moiety_names = getMoieties(moiety_fid)
  # Create the configuration file
  config_dct = {} 
  simple = SimpleSBML()
  simple.initialize(xml_fid)
  for molecule in simple.molecules:
    moiety_stoichiometries = findMoietyStoichiometries(
       molecule.name, moiety_names)
    config_dct[molecule.name] =  \
        [str(ms) for ms in moiety_stoichiometries]
  dct = {"moiety_structure": config_dct}
  yaml.dump(dct, config_fid)
  for fid in [xml_fid, moiety_fid, config_fid]:
    fid.close()

    
if __name__ == '__main__':
  def open_write(path):
    return open(path, "w")
  #
  parser = argparse.ArgumentParser(
      description='Creates configuration file for BioModels 147.')
  parser.add_argument('xml_file', 
      type=open, 
      help='SBML file from which species names are obtained.')
  parser.add_argument('moiety_file', 
      type=open, 
      help='YAML file the specifies moieties for the XML model.')
  parser.add_argument('--output', type=open_write,
      help="Generated configuration file. Defaults to config.yml.")
  args = parser.parse_args()
  if args.output is None:
    config_fid = open(DEFAULT_CONFIG_FILE, "w")
  else:
    config_fid = args.output
  main(args.xml_file, args.moiety_file, config_fid)
