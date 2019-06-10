"""Provides comparisons of moieties."""

from SBMLLint.common import constants as cn
from SBMLLint.common import config
from SBMLLint.common.reaction import Reaction
from SBMLLint.common.molecule import Molecule, MoleculeStoichiometry
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common.moiety  \
    import Moiety, MoietyStoichiometry

import collections
import pandas as pd
import numpy as np

NULL_STR = ''
INDENT = "  "


MoietyComparatorResult = collections.namedtuple('MoietyComparatorResult',
    'num_reactions num_imbalances report')


class MoietyComparator(object):
  """Analysis of moieties of molecules."""

  def __init__(self, mol_stoichs1, mol_stoichs2, 
      names=["reactants", "products"],
      implicits=None):
    """
    :param set-MoleculeStoichiometry mol_stoichs1:
    :param set-MoleculeStoichiometry mol_stoichs2:
    :param list-str names: names to refer to the two sets
    :param list-moieties implicits: implicit moieties
    """
    def checkType(objs):
      trues = [isinstance(o, MoleculeStoichiometry)
          for o in objs]
      if all(trues):
        return
      raise ValueError("Argument must be collection of type %s"
          % str(MoleculeStoichiometry))
    checkType(mol_stoichs1)
    checkType(mol_stoichs2)
    self.molecule_stoichiometry_collections = [
        mol_stoichs1, 
        mol_stoichs2,
        ]
    self.names = names
    if implicits is None:
      self._implicits = []
    else:
      self._implicits = implicits

  def _makeDFS(self):
    dfs = []
    for collection in self.molecule_stoichiometry_collections:
      if len(collection) > 0:
        df = MoleculeStoichiometry.countMoietysInCollection(collection)
        dfs.append(df)
    if len(dfs) == 0:
      dfs = [pd.DataFrame(), pd.DataFrame()]
    elif len(dfs) == 1:
      col = dfs[0].columns.tolist()[0]
      df = dfs[0].copy()
      df[col] = 0     
      dfs.append(df)
    return dfs

  def isSame(self):
    """
    Determines if the two molecules have the same type and count
    of moieties.
    :return bool:
    """
    dfs = self._makeDFS()
    return dfs[0].equals(dfs[1])

  def difference(self):
    """
    Calculates the moieties present in one set and not in the other.
    :return pd.DataFrame: row index is moiety,
        column Value: counts of moieties present in 0 
        and not in second. Negative values are present in 1 and
        not in 0.
    Handles boundary reactions, when one set of moieties
    is all zeroes.
    """
    def addDFIndex(df, index):
      # Adds missing indices to df
      missing = set(index).difference(df.index)
      for item in missing:
        df.loc[item] = 0
    def isZeroColumn(df):
      col = df.columns[0]
      return df[col].sum() == 0
    #
    dfs = self._makeDFS()
    addDFIndex(dfs[0], dfs[1].index)
    addDFIndex(dfs[1], dfs[0].index)
    drops = set(self._implicits).intersection(dfs[0].index)
    df0 = dfs[0].drop(drops)
    df1 = dfs[1].drop(drops)
    df = df0 - df1
    # Handle boundaries
    config_dict = config.getConfiguration()
    if not config_dict[cn.CFG_PROCESS_BOUNDARY_REACTIONS]:
      if isZeroColumn(df0) or isZeroColumn(df1):
        df[df.columns[0]] = 0
    #
    return df

  def reportDifference(self):
    """
    Reports a difference in moieties between the two sets of molecules.
    :return str: report. null if no difference.
    """
    if self.isSame():
      return NULL_STR
    df = self.difference()
    #
    def appendNewline(stg):
      if len(stg) > 0:
        return "%s\n" % stg
      else:
        return stg
    #
    def buildStg(name, sign, indent=INDENT):
      """
      Constructs the report for the name with the given
      sign of values in df.
      :param str name:
      :param int sign: 1 or -1
      :param str indent: indentation on each line
      :return str:
      """
      stg = NULL_STR
      for idx, row in df.iterrows():
        if sign*row[cn.VALUE] > 0:
          stg = "%s%s%s: %2.2f\n" % (
              stg, indent, idx, sign*row[cn.VALUE])
      if not stg == NULL_STR:
        stg = "Excess moieties in %s\n%s" % (name, stg)
      return stg
    #
    stg1 = appendNewline(buildStg(self.names[0], 1))
    stg2 = appendNewline(buildStg(self.names[1], -1))
    return "%s%s" % (stg1, stg2)

  @classmethod
  def analyzeReactions(cls, model_reference, implicits=None):
    """
    Analyzes all reactions to detect moiety imbalances.
    :param libsbml.Model or SimpleSBML model:
    :return MoietyComparatorResult:
    If model_reference is SimpleSBML, it must have been initialized.
    """
    if isinstance(model_reference, SimpleSBML):
      simple = model_reference
    else:
      simple = SimpleSBML()
      simple.initialize(model_reference)
    num_imbalances = 0
    report = NULL_STR
    for reaction in simple.reactions:
      comparator = cls(reaction.reactants, reaction.products,
          implicits=implicits)
      stg = comparator.reportDifference()
      if len(stg) > 0:
        num_imbalances += 1
        report = "%s\n***%s\n%s" % (
            report, 
            reaction.getId(is_include_kinetics=False),
            stg
            )
    num_reactions = len(simple.reactions)
    report = "\n%d of %d reactions have imbalances.\n%s" % (
        num_imbalances, num_reactions, report)
    result = MoietyComparatorResult(
        num_reactions=num_reactions,
        num_imbalances=num_imbalances,
        report=report)
    return result
