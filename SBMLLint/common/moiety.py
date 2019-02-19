"""
Implements a moiety and a moiety stoichiometry. A moiety
is a functional group in a molecule.
A moiety stoichiometry is a combination of a name and an integer
that indicates the repetitions of the moiety.
A moiety is separated from its stoichiometry by MOIETY_SEPARATOR.

Examples of 

   MOIETY STOICHIOMETRY STRING  MOIETY  STOICHIOMETRY
    A_1                          A      1
    A                            A      1
    A_P                          INVALID
"""

from SBMLLint.common import constants as cn
from SBMLLint.common import util


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


class MoietyStoichiometry(object):
  """A Moiety with its replication count."""

  def __init__(self, moiety, stoichiometry):
    if isinstance(moiety, Moiety):
      self.moiety = moiety
    else:
      self.moiety = Moiety(str(moiety))
    self.stoichiometry = stoichiometry
    self.name = "%s%s%d" % (self.moiety.name,
        cn.MOIETY_SEPARATOR, stoichiometry)

  def __repr__(self):
    return self.name

  def __lt__(self, other):
    return str(self) < str(other)

  def isEqual(self, other):
    return self.moiety.isEqual(other.moiety) and  \
        (self.stoichiometry == other.stoichiometry)

  @classmethod
  def getMoietys(cls, moiety_stoichiometrys):
    """
    Extract moieties from MoietyStoichiometrys
    """
    moietys = util.uniqueify([m_s.moiety 
        for m_s in moiety_stoichiometrys])
    return moietys

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
    return cls(Moiety(name), stoich)
