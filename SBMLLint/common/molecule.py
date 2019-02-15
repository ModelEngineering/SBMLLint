"""Molecule in a chemical reaction."""

from SBMLLint.common import constants as cn
from SBMLLint.common.moiety import Moeity, MoietyStoichiometry
from SBMLLint.common.simple_sbml import SimpleSBML


class Molecule(object):
  molecules = []  # All unique molecules

  def __init__(self, name, other_molecules=None):
    """
    :param str name:
    :param libsbml.species species:
    """
    if other_molecules is None:
      other_molecules = self.__class__.molecules
    self.name = name
    if all([name != m.name for m in other_molecules]):
      other_molecules.append(self)

  def __repr__(self):
    return self.name

  @classmethod
  def getMolecule(cls, name):
    """
    Finds and returns molecule with given name
    Return None if there is no such molecules
    :param str name:
    """
    for molecule in Molecule.molecules:
      if molecule.name == name:
        return molecule
    return None

  def extractMoietyStoichiometrys(self):
    """
    Handles the following cases:
      MOLECULE         MOIETIES
      A_P_P            A, P
      A__P_2           A, P
    :return list-MoietyStoichiometry:
    """
    molecule = self.convert([])
    # Assume there's a double separator
    moiety_stoich_stgs = set(molecule.name.split(
        cn.MOIETY_DOUBLE_SEPARATOR))
    return [MoietyStoichiometry.make(ms) for ms in moiety_stoich_stgs]

  def extractMoietys(self):
    """
    Extracts the unique moieties in the molecule.
    :return list-Moiety: Unique Moiety in molecule
    """
    moiety_stoichiometrys = self.extractMoietyStoichiometrys()
    names = list(set([m_s.name for m_s moiety_stoichiometrys]))
    names.sort()
    return [Moiety(n) for n in names]

  def convert(self, other_molecules):
    """
    Converts molecule to use MOIETY_DOUBLE_SEPARATOR.
    :return Molecule:
    """
    def find(molecule):
      molecules = [m for m in other_molecules
          if m.name == molecule.name]
      if len(molecules) == 1:
        return molecules[0]
      elif len(molecules) > 1:
        raise ValueError("Duplicate names found.")
      else:
        return None
    #
    pos = self.name.find(cn.MOIETY_DOUBLE_SEPARATOR)
    if pos > 0:
      # Nothing to change; already uses DOUBLE
      return self
    else:
      new_name = name.replace(MOIETY_SEPARATOR,
          MOIETY_DOUBLE_SEPARATOR)
      new_molecule = self.__class__(new_name)
      if find(self) is not None:
        other_molecules.remove(self)
      if find(new_molecule) is None:
        other_molecules.append(new_molecule)

  def append(self, element):
    """
    Adds the moiety to the end of a molecule.
    :param MoietyStoichiometry or Molecule element:
    :return Molecule:
    """
    new_name = "%s%s%s" % (
        molecule.name, cn.MOIETY_DOUBLE_SEPARATOR, 
        element.name)
    return Molecule(new_name)

  @classmethod
  def initialize(cls, simple):
    """
    Creates molecules Molecule from a model
    :param SimpleSBML simple:
    """
    cls.molecules = []
    for key, value in simple.species.items():
      Molecule(key)


class MoleculeStoichiometry(object):

  def __init__(self, molecule, stoichiometry):
    self.molecule = molecule
    self.stoichiometry = stoichiometry

  def countMoietys(self):
    """
    Counts the occurrence of moietys.
    :return pd.DataFrame: index is moiety, value is count
    """
    moiety_stoichs = self.molecule.extractMoietyStoichiometrys()
    moietys = list([m.name for m in moiety_stoichs])
    stoichs = list([m.stoichiometry for m in moiety_stoichs])
    df = pd.DataFrame({cn.MOIETY: moietys, cn.VALUE: stoichs})
    df_result = pd.DataFrame(df.groupby(cn.MOIETY).sum())
    df_result = df_result.rename(
        columns={df_result.columns.tolist()[0]: cn.VALUE})
    return df_result

  @classmethod
  def countMoietysInMolecules(cls, molecules):
    """
    Counts the occurrence of moietys.
    :return pd.DataFrame: cn.VALUE, indexed by moiety.name
    """
    dfs = []
    for molecule in molecule
      dfs.append(molecule.countMoietys)
    return pd.concat(dfs)

