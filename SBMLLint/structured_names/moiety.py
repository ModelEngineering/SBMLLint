"""
Implements a moiety, part of a molecule, and a moiety stoichiometry.
A moiety stoichiometry is a combination of a name and number
that indicates how many times the moiety occurs. If the number
is present, moieties in a molecule are 
separated by a MOIETY_DOUBLE_SEPARATOR;
otherwise, moieties can be separated by a MOIETY_SEPARATOR.

Examples:

   MOLECULE       LIST OF MOIETY, NUMBER
    A_P            (A, 1), (P, 1)
    A__P           (A, 1), (P, 1)
    A_1__P_1       (A, 1), (P, 1)
    A_1__P         (A, 1), (P, 1)
    A_1.0__P_1     (A, 1.0), (P, 1)
    A_1_P          INVALID
"""

from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common.simple_sbml import SimpleSBML

import pandas as pd
import numpy as np

NULL_STR = ''

############## INERNAL FUNCTIONS ##################

def _extractFromMoietyStoichiometryString(moiety_stoich_stg):
  """
  Extracts the name and count from a string for a
  moiety stoichiometry. Examples are:
     STRING     RESULT
      P_2        P, 2
      A          A, 1
  """
  result = []
  pos = moiety_stoich_stg.find(cn.MOIETY_SEPARATOR)
  if pos < 0:
    name = moiety_stoich_stg
    stoich_stg = 1
  elif pos == 0:
    raise ValueError(
        "Invalid format for moiety stoichiometry string: %s"
        % moiety_stoich_stg)
  else:
    name = moiety_stoich_stg[0:pos]
    stoich_stg = moiety_stoich_stg[pos+1:]
  try:
    stoich = float(stoich_stg)
  except ValueError:
    raise ValueError(
        "Invalid number in moiety stoichiometry string: %s"
        % moiety_stoich_stg)
  return name, stoich


############## CLASSES ##################
class Moiety(object):

  def __init__(self, name, other_moietys=[]):
    """
    :param str name:
    :param list-Moiety other_moeitys:
    Ensures unique names within other_moietys
    """
    self.name = name
    if all([name != m.name for m in other_moietys]):
      other_moietys.append(self)

  def __repr__(self):
    return self.name

  def __lt__(self, other):
    """
    Enables sorting a list of Moiety
    """
    return self.name < other.name

  def isEqual(self, other):
    return self.name == other.name

  @classmethod
  def extract(cls, molecule):
    """
    Creates moieties present in molecule.
    Handles the following cases:
      MOLECULE         MOIETIES
      A_P_P            A, P
      A__P_2           A, P
    :param Molecule molecule:
    :return list-Moiety: Unique moieties found
    """
    # Assume there's a double separator
    moiety_stoich_stgs = set(molecule.name.split(
        cn.MOIETY_DOUBLE_SEPARATOR))
    # See if double separator present
    if len(moiety_stoich_stgs) <= 1:
      moiety_stoich_stgs = set(molecule.name.split(cn.MOIETY_SEPARATOR))
    names = list(set([_extractFromMoietyStoichiometryString(n)[0] 
        for n in moiety_stoich_stgs]))
    names.sort()
    return [Moiety(n) for n in names]

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
    """
    :return bool: True if moiety name in molecule
    """
    moietys = molecule.name.split(cn.MOIETY_SEPARATOR)
    return self.name in moietys


class MoietyStoichiometry(object):
  """A Moiety with its count."""

  def __init__(self, moiety, stoichiometry):
    self.moiety = moiety
    self.stoichiometry = stoichiometry

  def __repr__(self):
    return "%s: %2.2f" % (str(self.moiety), self.stoichiometry)

  @classmethod
  def extract(cls, mole_stoich):
    """
    Creates moieties and their stoichiometries.
    :param MoleculeStoichiometry mole_stoich:
    :return list-MoietyStoichiometry:
    """
    names = set(mole_stoich.molecule.name.split(cn.MOIETY_SEPARATOR))
    result = []
    for name in names:
      num_occurrences = mole_stoich.molecule.name.count(name)
      stoich = mole_stoich.stoichiometry*num_occurrences
      new_moiety_stoich = MoietyStoichiometry(name, stoich)
      result.append(new_moiety_stoich)
    return result

  @classmethod
  def countMoietys(cls, mole_stoichs):
    """
    Counts the occurrence of moietys in a 
    list of MoleculeStoichiometry.
    :param list-MoleculeStoichiometry mole_stoichs:
    :return pd.DataFrame: index is moiety, value is count
    """
    moiety_stoichs = []
    for mole_stoich in mole_stoichs:
      moiety_stoichs.extend(cls.extract(mole_stoich))
    moietys = list(set([m.molecule.name for m in mole_stoichs]))
    counts = list(set([m.stoichiometry for m in mole_stoichs]))
    df = pd.DataFrame({cn.MOIETY: moietys, cn.VALUE: counts})
    df_result = pd.DataFrame(df.groupby(cn.MOIETY).sum())
    df_result = df_result.rename(
        columns={df_result.columns.tolist()[0]: cn.VALUE})
    return df_result
