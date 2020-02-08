#!/usr/bin/env python
"""
Runs moiety analysis for a local XML file.
Usage: moiety_analysis <filepath>
"""

from SBMLLint.common import constants as cn
from SBMLLint.tools import sbmllint

import argparse

def main():
  parser = argparse.ArgumentParser(description='SBML XML file.')
  parser.add_argument('xml_file', type=open, help='SBML file')
  parser.add_argument('--config', type=open,
      help="SBMLLint configuration file")
  args = parser.parse_args()
  sbmllint.lint(model_reference=args.xml_file,
      mass_balance_check=cn.MOIETY_ANALYSIS,
      config_fid=args.config)


if __name__ == '__main__':
  main()
