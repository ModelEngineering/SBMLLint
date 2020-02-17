#!/usr/bin/env python
"""
Runs the GAMES algorithm for a local XML file.
Usage: games <filepath>
"""

from SBMLLint.common import util
from SBMLLint.tools import sbmllint

import argparse

def main():
  parser = argparse.ArgumentParser(description='SBML XML file.')
  parser.add_argument('xml_file', type=open, help='SBML file')
  parser.add_argument('--config', type=open,
      help="SBMLLint configuration file")
  args = parser.parse_args()
  for fid in util.getNextFid(args.xml_file):
    try:
      sbmllint.lint(model_reference=fid,
          config_fid=args.config, mass_balance_check=sbmllint.GAMES)
    except ValueError:
      print ("  *** Bad SBML file.")


if __name__ == '__main__':
  main()
