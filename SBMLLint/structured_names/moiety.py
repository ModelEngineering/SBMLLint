"""Implements a moiety, part of a molecule."""

from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common.simple_sbml import SimpleSBML

import pandas as pd
import numpy as np


class Moiety(object):
  moietys = []  # All distinct moietys

  def __init__(self, name):
    self.name = name

  @classmethod
  def extract(cls, molecule):
    """
    Creates the moietys from a molecule.
    Maintains list of moietys.
    :param Molecule molecule:
    :return list of moietys:
    """
    names = molecule.name.split(cn.MOIETY_SEPARATOR)
    result = [cls(n) for n in names]
    cls.moietys.append(list(set(result)))
    return result

  def appendToMolecule(self, molecule):
    """
    Adds the moiety to the end of a molecule.
    :param Molecule molecule:
    :return Molecule:
    """
    new_name = "%s%s%s" % (
        molecule.name, cn.MOIETY_SEPARATOR, self.name)
    return Molecule(new_name)

  def isInMolecule(self, molecule):
    moietys = molecule.name.split(cn.MOIETY_SEPARATOR)
    return self.name in moietys

  @classmethod
  def countMoietys(cls, molecules):
    """
    Counts the occurrence of moietys in a list of molecules.
    :param list-Molecule molecules:
    :return pd.DataFrame: index is moiety, value is count
    """
    moietys = []
    for molecule in molecules:
      moietys.extend([m.name for m in cls.extract(molecule)])
    df = pd.DataFrame({cn.VALUE: moietys})
    return pd.DataFrame(df.groupby(cn.VALUE).size())


class MoietyComparator(object):
  """Analysis of moieties of molecules."""

  def __init__(self, molecules1, molecules2):
    """
    :param set-Molecule molecules1:
    :param set-Molecule molecules2:
    """
    self.list_of_molecules = [molecules1, molecules2]

  def isSame(self):
    """
    Determines if the two molecules have the same type and count
    of moieties.
    :return bool:
    """
    df0 = Moiety.countMoietys(self.list_of_molecules[0])
    df1 = Moiety.countMoietys(self.list_of_molecules[1])
    return df0.equals(df1)

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
    df0 = Moiety.countMoietys(self.list_of_molecules[0])
    df1 = Moiety.countMoietys(self.list_of_molecules[1])
    addDFIndex(df0, df1.index)
    addDFIndex(df1, df0.index)
    return df0 - df1
