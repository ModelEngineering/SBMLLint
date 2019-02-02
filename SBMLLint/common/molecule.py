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

  @classmethod
  def addMolecule(cls, molecule):
    if any([m.name == molecule.name for m in cls.molecules]):
      pass
    else:
      cls.molecules.append(molecule)

  @classmethod
  def initialize(cls, simple):
    """
    Creates molecules Molecule from a model
    :param SimpleSBML simple:
    """
    cls.molecules = []
    for key, value in simple.species.items():
      Molecule(key, species=value)

