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
    A_1_P          INVALID
"""

from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common.simple_sbml import SimpleSBML

import pandas as pd
import numpy as np

NULL_STR = ''



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

  def isInMolecule(self, molecule):
    """
    :return bool: True if moiety name in molecule
    """
    moietys = molecule.name.split(cn.MOIETY_SEPARATOR)
    return self.name in moietys


class MoietyStoichiometry(object):
  """A Moiety with its replication count."""

  def __init__(self, moiety, stoichiometry):
    self.moiety = moiety
    self.stoichiometry = stoichiometry
    self.name = "%s%s%d" % (self.moiety.name,
        cn.MOIETY_SEPARATOR, stoichiometry)

  def __repr__(self):
    return self.name

  @classmethod
  def make(cls, moiety_stoich_stg):
    """
    Makes a MoietyStoichiometry from a string.
    Examples of strings are: "P_2", "A"
    :return MoietyStoichiometry:
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
      stoich = int(stoich_stg)
    except ValueError:
      raise ValueError(
          "Invalid number in moiety stoichiometry string: %s"
          % moiety_stoich_stg)
    return cls(name, stoich)
