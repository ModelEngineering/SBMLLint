"""Molecule in a chemical reaction."""

from SBMLLint.common import constants as cn
from SBMLLint.common.simple_sbml import SimpleSBML


class Molecule(object):
  all = []  # All molecules

  def __init__(self, name, species=None):
    """
    :param str name:
    :param libsbml.species species:
    """
    self.name = name
    self._species = species

  @classmethod
  def initialize(cls, simple):
    """
    Creates all Molecule from a model
    :param SimpleSBML simple:
    """
    cls.all = []
    for key, value in simple.species.items():
      cls.all.append(Molecule(key, species=value))

