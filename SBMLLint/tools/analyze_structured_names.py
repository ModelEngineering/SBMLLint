"""A Tool that creates data to analyze the use of structured names."""

from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common import simple_sbml
from SBMLLint.common import util
from SBMLLint.common.reaction import Reaction
from SBMLLint.tools import sbmllint
from SBMLLint.tools import print_reactions

import os
import pandas as pd

# Columns in output file
FILENAME = "filename"
HAS_STRUCTURE = "has_structure"
NUM_REACTIONS = "num_reactions"
NUM_BOUNDARY_REACTIONS = "num_boundary_reactions"
NUM_BAD = "num_imbalance_reactions"
NUM_BALANCED_REACTIONS = "num_balanced_reactions"
FRC_BALANCED = "frc_balanced"
# File paths
OUTPUT_FILE = "analyze_structured_names.csv"
DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_PATH = os.path.join(DIR, OUTPUT_FILE)
# Miscellaneous
EXCLUDE_LIST = ["node", "x", "species"]

def isStructuredName(name):
  """
  Detects if the name is a candidate for a structured name.
  The term "element" refers to substrings separated by "_".
  Not a structured name if any of the following are True:
    No "_" present
    Last element is a number
    First element is in EXCLUDE_LIST (case independent)
  """
  moietys = name.split(cn.MOIETY_SEPARATOR)
  if len(moietys) == 1:
    return False
  # See if last element has a numeric suffix
  try:
    _ = float(moietys[-1])
    return False
  except:
      pass
  if any([moietys[0].lower() in s for s in EXCLUDE_LIST]):
    return False
  return True

def calcStats(initial=1, final=50, out_path=OUTPUT_PATH, 
    report_interval=50, report_progress=True, min_frc=-1):
  """
  Calculates statistics for structured names.
  :param int initial: Index of first model to process
  :param int final: Index of final model to process
  :param str out_path: Path to the output CSV file
  :param int report_interval: Number of files processed before
      a report is written
  :param bool report_progress: report file being processed
  :param float min_frc: Filter to select only those models
      that have at least the specified fraction of reactions
      balanced according to structured_names
  """
  def writeDF(dfs):
    df_count = pd.concat(dfs)
    df_count[NUM_BALANCED_REACTIONS] = df_count[NUM_REACTIONS] - df_count[NUM_BAD]
    df_count[FRC_BALANCED] = 1.0*df_count[NUM_BALANCED_REACTIONS] / (
        df_count[NUM_REACTIONS] - df_count[NUM_BOUNDARY_REACTIONS])
    df = df_count[df_count["frc_balanced"] > min_frc]
    df = df.sort_values("frc_balanced")
    df.to_csv(out_path, index=False)
  #
  dfs = []
  sbmliter = simple_sbml.modelIterator(initial=initial, final=final)
  for item in sbmliter:
    if report_progress:
      print("*Processing file %s, number %d"
           % (item.filename, item.number))
    simple = simple_sbml.SimpleSBML()
    try:
      simple.initialize(item.model)
    except:
      print("  Error in model number %d." % item.number)
      continue
    row = {FILENAME: [item.filename], 
           HAS_STRUCTURE: [False], 
           NUM_BOUNDARY_REACTIONS: [0],
           NUM_REACTIONS: [0],
           NUM_BAD: [None],
           NUM_BAD: [0],
           }
    for reaction in simple.reactions:
      if (len(reaction.reactants) == 0) or (len(reaction.products) == 0):
          row[NUM_BOUNDARY_REACTIONS] = [row[NUM_BOUNDARY_REACTIONS][0] + 1]
      molecules = util.uniqueify([m.molecule 
          for m in set(reaction.reactants).union(reaction.products)])
      if any([isStructuredName(m.name) for m in molecules]):
          row[HAS_STRUCTURE] = [True]
    try:
      mcr = sbmllint.lint(item.model, is_report=False)
      row[NUM_REACTIONS] = [mcr.num_reactions]
      row[NUM_BAD] = [mcr.num_imbalances]
    except:
      row[NUM_REACTIONS] = [None]
      row[NUM_BAD] = [None]
    dfs.append(pd.DataFrame(row))
    if item.number % report_interval == 0:
      writeDF(dfs)
  writeDF(dfs)


if __name__ == '__main__':
  calcStats(initial=1, final=1000, report_interval=1)
