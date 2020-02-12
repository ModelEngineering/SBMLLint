"""A Tool that creates data to analyze the use of names
   that are structured to expose moieties.
"""

from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common import simple_sbml
from SBMLLint.common import util
from SBMLLint.common.reaction import Reaction
from SBMLLint.tools import sbmllint
from SBMLLint.tools import print_reactions

import os
import numpy as np
import pandas as pd

# File paths
OUTPUT_FILE = "analyze_moiety_analysis.csv"
DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_PATH = os.path.join(DIR, OUTPUT_FILE)
# Miscellaneous
EXCLUDE_PREFIX = ["node", "x", "species"]
EXCLUDE_SUFFIX = ["prime", "mrna", "prot"]

def isStructuredName(name):
  """
  Detects if the name is a candidate for a moiety structured name.
  The term "element" refers to substrings separated by "_".
  Not a structured name if any of the following are True:
    No "_" present
    Last element is a number
    First element is in EXCLUDE_PREFIX (case independent)
    Last element is in EXCLUDE_SUFFIX (case independent)
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
  if any([moietys[0].lower() in s for s in EXCLUDE_PREFIX]):
    return False
  if any([moietys[-1].lower() in s for s in EXCLUDE_SUFFIX]):
    return False
  return True

def calcStats(initial=0, final=50, out_path=OUTPUT_PATH, 
    report_interval=50, report_progress=True, min_frc=-1,
    data_dir=cn.BIOMODELS_DIR):
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
      balanced according to moiety_analysis
  """
  def writeDF(dfs):
    df_count = pd.concat(dfs)
    df_count[cn.NUM_BALANCED_REACTIONS] =  \
        df_count[cn.TOTAL_REACTIONS]  \
        - df_count[cn.NUM_IMBALANCED_REACTIONS]
    denom =  (df_count[cn.TOTAL_REACTIONS] 
        - df_count[cn.NUM_BOUNDARY_REACTIONS])
    denom = [np.nan if np.isclose(v, 0) else v for v in denom]
    df_count[cn.FRAC_BALANCED_REACTIONS] =  \
        1.0*df_count[cn.NUM_BALANCED_REACTIONS] / denom
    df_count[cn.FRAC_BOUNDARY_REACTIONS] =  \
        1.0*df_count[cn.NUM_BOUNDARY_REACTIONS] / (
        df_count[cn.TOTAL_REACTIONS])
    if min_frc < 0:
      df = df_count
    else:
      df = df_count[df_count[cn.FRAC_BALANCED_REACTIONS] > min_frc]
    df = df.sort_values(cn.FRAC_BALANCED_REACTIONS)
    df.to_csv(out_path, index=False)
  #
  dfs = []
  sbmliter = simple_sbml.modelIterator(initial=initial, final=final,
      data_dir=data_dir)
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
    row = {cn.FILENAME: [item.filename], 
           cn.IS_STRUCTURED: [False], 
           cn.NUM_BOUNDARY_REACTIONS: [0],
           cn.TOTAL_REACTIONS: [0],
           cn.NUM_IMBALANCED_REACTIONS: [0],
           }
    for reaction in simple.reactions:
      if (len(reaction.reactants) == 0) or (len(reaction.products) == 0):
          row[cn.NUM_BOUNDARY_REACTIONS] =  \
              [row[cn.NUM_BOUNDARY_REACTIONS][0] + 1]
      molecules = util.uniqueify([m.molecule 
          for m in set(reaction.reactants).union(reaction.products)])
      if any([isStructuredName(m.name) for m in molecules]):
          row[cn.IS_STRUCTURED] = [True]
    try:
      mcr = sbmllint.lint(model_reference=item.model, is_report=False)
      row[cn.TOTAL_REACTIONS] = [mcr.num_reactions if mcr.num_reactions > 0 else np.nan]
      row[cn.NUM_IMBALANCED_REACTIONS] = [mcr.num_imbalances]
    except:
      row[cn.TOTAL_REACTIONS] = [None]
      row[cn.NUM_IMBALANCED_REACTIONS] = [0]
    dfs.append(pd.DataFrame(row))
    if item.number % report_interval == 0:
      writeDF(dfs)
  writeDF(dfs)


if __name__ == '__main__':
  calcStats(initial=0, final=1000, report_interval=1,
      data_dir=cn.BIGG_DIR)
