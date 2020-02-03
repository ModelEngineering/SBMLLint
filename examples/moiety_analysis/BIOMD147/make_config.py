import argparse
import yaml
from SBMLLint.common.simple_sbml import SimpleSBML

config_start = """
# Configuration for SBMLLint
# Implicit molecules are molecules that may be implicit
# in either the reactants or products.
implicits:
- cytoplasm
- nucleus
- IkBa_mRNA
- IkBe_mRNA
- IkBa
- IkBb
- IkBe

# Indicates if boundary reactions are checked
process_boundary_reactions: False

# Explicit declaration of moiety structures
"""

# Molecules contain one or more of the names in MOIETIES,
# and moiety stoichiometries are always 1.
MOIETIES = ["IkBa", "IkBb", "IkBe", "IKK", 
             "NFkB", "nucleus", "cytoplasm"]
DEFAULT_CONFIG_FILE = "config.yml"

def main(xml_fid, config_file):
  """
  Creates the YAML configuration file.
  :param IOWrappter xml_fid: xml file
  :param str config_file: file created
  """
  config_dct = {} 
  simple = SimpleSBML()
  simple.initialize(xml_fid)
  for molecule in simple.molecules:
    config_dct[molecule.name] =  \
        ["%s, 1" % m for m in MOIETIES if m in molecule.name]
  dct = {"moiety_structure": config_dct}
  with open(config_file, "w") as fd:
    fd.writelines(config_start)
    yaml.dump(dct, fd)

    
if __name__ == '__main__':
  parser = argparse.ArgumentParser(
      description='Creates configuration file for BioModels 147.')
  parser.add_argument('xml_file', 
      type=open, 
      help='SBML file from which species names are obtained.')
  parser.add_argument('--output', type=str,
      help="Generated configuration file. Defaults to config.yml.")
  args = parser.parse_args()
  if args.output is None:
    config_file = DEFAULT_CONFIG_FILE
  else:
    config_file = args.output
  main(args.xml_file, config_file)
