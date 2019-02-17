"""Provides comparisons of moieties."""

from SBMLLint.common import constants as cn
from SBMLLint.common.reaction import Reaction
from SBMLLint.common.molecule import Molecule, MoleculeStoichiometry
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common.moiety  \
    import Moiety, MoietyStoichiometry

import pandas as pd
import numpy as np

NULL_STR = ''
INDENT = "  "


class MoietyComparator(object):
  """Analysis of moieties of molecules."""

  def __init__(self, mol_stoichs1, mol_stoichs2, 
      names=["reactants", "products"]):
    """
    :param set-MoleculeStoichiometry mol_stoichs1:
    :param set-MoleculeStoichiometry mol_stoichs2:
    :param list-str names: names to refer to the two sets
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

  def _makeDFS(self):
    dfs = []
    for collection in self.molecule_stoichiometry_collections:
      df = MoleculeStoichiometry.countMoietysInCollection(collection)
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
    """
    def addDFIndex(df, index):
      # Adds missing indices to df
      missing = set(index).difference(df.index)
      for item in missing:
        df.loc[item] = 0
    #
    dfs = self._makeDFS()
    addDFIndex(dfs[0], dfs[1].index)
    addDFIndex(dfs[1], dfs[0].index)
    return dfs[0] - dfs[1]

  def reportDifference(self):
    """
    Reports a difference in moieties between the two sets of molecules.
    :return str: report. null if no difference.
    """
    if self.isSame():
      return NULL_STR
    df = self.difference()
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
  def initialize(cls, simple):
    Reaction.initialize(simple)

  @classmethod
  def analyzeReactions(cls, model):
    """
    Analyzes all reactions to detect moiety imbalances.
    :param libsbml.Model or SimpleSBML model:
    :return int, str: number of reactions with imbalances,
        report
    """
    if isinstance(model, SimpleSBML):
      simple = model
    else:
      simple = SimpleSBML(model)
    cls.initialize(simple)
    num = 0
    report = NULL_STR
    for reaction in Reaction.reactions:
      comparator = cls(reaction.reactants, reaction.products)
      stg = comparator.reportDifference()
      if len(stg) > 0:
        num += 1
        report = "%s\n***%s\n%s" % (
            report, 
            reaction.getId(is_include_kinetics=False),
            stg
            )
    report = "\n%d reactions have imbalances.\n%s" % (num, report)
    return num, report
