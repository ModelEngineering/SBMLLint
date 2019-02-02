"""Implements a moiety, part of a molecule."""

from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common.simple_sbml import SimpleSBML

import pandas as pd
import numpy as np


class Moiety(object):
  moietys = []  # All distinct moietys

  def __init__(self, name)
    self._name = name

  @classmethod
  def make(cls, molecule):
    """
    Creates the moietys from a molecule.
    Maintains list of moietys.
    :param Molecule molecule:
    :return list of moietys:
    """
    names = molecule.name.split(cn.MOIETY_SEPARATOR)
    result = [cls(n) for n in names]
    cls.moietys.append(list(set(result)))

  @classmethod
  def addMoietyToMolecule(cls, molecule, moiety):
    """
    Adds the moiety to a molecule.
    :param Molecule molecule:
    :param Moiety moeity:
    :return Molecule:
    """
    new_name = "%s%s%s" % (
        molecule.name, cn.MOIETY_SEPARATOR, moiety.name)
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
      moietys.extend(cls.make(molecule))
    df = pd.DataFrame(cn.VALUE: moietys)
    return pd.DataFrame(df.groupby(cn.VALUE).size())
    
