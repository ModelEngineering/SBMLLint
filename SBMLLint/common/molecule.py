"""Molecule in a chemical reaction."""

from SBMLLint.common import constants as cn
from SBMLLint.common.simple_sbml import SimpleSBML


class Molecule(object):
  molecules = []  # All unique molecules

  def __init__(self, name, species=None):
    """
    :param str name:
    :param libsbml.species species:
    """
    self.name = name
    self._species = species
    self.__class__.addMolecule(self)

  def __repr__(self):
    return self.name

  @classmethod
  def addMolecule(cls, molecule):
    if any([m.name == molecule.name for m in cls.molecules]):
      pass
    else:
      cls.molecules.append(molecule)


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


  @classmethod
  def initialize(cls, simple):
    """
    Creates molecules Molecule from a model
    :param SimpleSBML simple:
    """
    cls.molecules = []
    for key, value in simple.species.items():
      Molecule(key, species=value)

