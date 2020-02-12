"""
Molecule in a chemical reaction and MoleculeStoichiometry
(a molecule with its stoichiometry).

A Molecule is structured as one or more MoietyStoichiometry
with a separator. MOIETY_DOUBLE_SEPARATOR is used if at least
one MoietyStoichiometry has a repetition count; otherwise,
either MOIETY_DOUBLE_SEPARATOR or MOIETY_SEPRATOR can be used.

MOLECULE           MOIETY, STOICHIOMETRY
 A                  (A, 1)
 A_P_P_P            (A, 1), (P, 3)
 A__P__P__P         (A, 1), (P, 3)
 A__P_3             (A, 1), (P, 3)
"""

from SBMLLint.common import constants as cn
from SBMLLint.common import config
from SBMLLint.common.moiety import Moiety, MoietyStoichiometry
from SBMLLint.common import util

import pandas as pd
import numpy as np


class Molecule(object):

  def __init__(self, name):
    """
    :param str name:
    """
    self.name = name
    self._moiety_stoichiometrys = None

  @property
  def moiety_stoichiometrys(self):
    done = False
    if self._moiety_stoichiometrys is None:
      config_dct = config.getConfiguration()
      if cn.CFG_MOIETY_STRUCTURE in config_dct:
        dct = config_dct[cn.CFG_MOIETY_STRUCTURE]
        if self.name in dct.keys():
          self._moiety_stoichiometrys =  \
              MoietyStoichiometry.makeFromDct(dct[self.name])
          done = True
    else:
      done = True
    if not done:
      new_name = self._reformat()
      stgs = new_name.split(cn.MOIETY_DOUBLE_SEPARATOR)
      result = [MoietyStoichiometry.make(ms) for ms in stgs]
      result.sort()
      self._moiety_stoichiometrys = result
    return self._moiety_stoichiometrys

  def __repr__(self):
    return self.name

  def __lt__(self, other):
    return self.name < other.name

  def isEqual(self, other):
    return self.name == other.name

  def getMoietys(self):
    """
    Extracts the unique moieties in the molecule.
    :return list-Moiety: Unique Moiety in molecule
    """
    names = list(set([m_s.moiety.name for m_s in self.moiety_stoichiometrys]))
    names.sort()
    return [Moiety(n) for n in names]

  def _reformat(self):
    """
    Reformats the molecule name to use MOIETY_DOUBLE_SEPARATOR.
    :return str:
    """
    new_name = None
    pos = self.name.find(cn.MOIETY_DOUBLE_SEPARATOR)
    if pos > 0:
      new_name = self.name
    else:
      # Check to see if there is a single moiety
      parts = self.name.split(cn.MOIETY_SEPARATOR)
      if len(parts) == 2:
        if util.isInt(parts[1]):
          new_name = self.name
      if new_name is None:
        new_name = self.name.replace(cn.MOIETY_SEPARATOR,
            cn.MOIETY_DOUBLE_SEPARATOR)
    return new_name

  def append(self, element):
    """
    Appends to the end of a molecule.
    :param MoietyStoichiometry or Molecule element:
    :return Molecule:
    """
    new_name = "%s%s%s" % (
        self.name, cn.MOIETY_DOUBLE_SEPARATOR, 
        element.name)
    return Molecule(new_name)

  def hasMoiety(self, moiety):
    """
    :param Moiety moiety:
    :return bool: True if moiety name in molecule
    """
    moietys = self.getMoietys()
    return any([moiety.isEqual(m) for m in moietys])


class MoleculeStoichiometry(object):

  def __init__(self, molecule, stoichiometry):
    if not isinstance(molecule, Molecule):
      raise ValueError("First argument must be a Molecule.")
    if not util.isFloat(stoichiometry):
      raise ValueError("Second argument must be a float.")
    self.molecule = molecule
    self.stoichiometry = stoichiometry

  def __repr__(self):
    return "%s * % 2.2f" % (str(self.molecule), self.stoichiometry)

  def __lt__(self, other):
    return str(self) < str(other)

  def countMoietys(self):
    """
    Counts the occurrence of moietys.
    :return pd.DataFrame: index is moiety, value is count
    """
    moiety_stoichs = self.molecule.moiety_stoichiometrys
    moietys = list([str(m.moiety) for m in moiety_stoichs])
    stoichs = list([m.stoichiometry for m in moiety_stoichs])
    df = pd.DataFrame({cn.MOIETY: moietys, cn.VALUE: stoichs})
    df_result = pd.DataFrame(df.groupby(cn.MOIETY).sum())
    df_result = df_result.rename(
        columns={df_result.columns.tolist()[0]: cn.VALUE})
    df_result[cn.VALUE] = self.stoichiometry*df_result[cn.VALUE]
    return df_result

  @classmethod
  def countMoietysInCollection(cls, molecule_stoichiometrys):
    """
    Counts the occurrence of moietys.
    :param list-MoleculeStoichiometry  molecule_stoichiometrys:
    :return pd.DataFrame: cn.VALUE, indexed by moiety.name
    """
    dfs = []
    for molecule_stoichiometry in molecule_stoichiometrys:
      dfs.append(molecule_stoichiometry.countMoietys())
    df = pd.concat(dfs)
    df = df.reset_index()
    col = df.columns[0]
    df_result = pd.DataFrame(df.groupby(col).sum())
    return df_result

  @classmethod
  def getMolecules(cls, molecule_stoichiometrys):
     """
     Obtainins the molecules from a collection of
     MoleculeStoichiometry.
     :param list-MoleculeStoichiometry molecule_stoichiometrys;
     :return list-Molecule: unique molecules
     """
     molecules = [m_s.molecule for m_s in molecule_stoichiometrys] 
     return util.uniqueify(molecules)
