"""Chemical Reaction."""

from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common.simple_sbml import SimpleSBML

import numpy as np


REACTION_SEPARATOR = "->"


class Reaction(object):
  reactions = []  # All reactions

  def __init__(self, libsbml_reaction):
    self.identifier = SimpleSBML.getReactionString(libsbml_reaction,
        is_include_kinetics=True)
    self.reactants = self._getMolecules(libsbml_reaction,
        SimpleSBML.getReactants)
    self.products = self._getMolecules(libsbml_reaction,
        SimpleSBML.getProducts)
    self.category = self._getCategory()
    if not any([self.isEqual(r) for r in Reaction.reactions]):      
      self.__class__.reactions.append(self)

  def getId(self, is_include_kinetics=True):
    if is_include_kinetics:
      return self.identifier
    else:
      pos = self.identifier.find(cn.KINETICS_SEPARATOR)
      if pos < 0:
        return self.identifier
      else:
        return self.identifier[:pos]

  def __repr__(self):
    return self.getId()

  def _getMolecules(self, libsbml_reaction, func):
    """
    Constructs molecules for the species returned by function.
    :param Function func: gets libsbml.SpeciesReference
    :return list-MoleculeStoichiometry:
    """
    species = func(libsbml_reaction)
    molecules = []
    for spc in species:
      molecule = Molecule.getMolecule(spc.species)
      if molecule is None:
          molecule = Molecule(spc.species)
      molecules.append(cn.MoleculeStoichiometry(
          molecule = molecule,
          stoichiometry = spc.getStoichiometry())
          )
    return molecules

  def _getCategory(self):
    """
    :return str: reaction category
    """
    num_reactants = len(self.reactants)
    num_products = len(self.products)
    for reaction_category in cn.REACTION_CATEGORIES:
      if reaction_category.predicate(num_reactants, num_products):
        return reaction_category.category
    raise ValueError("Reaction category not found.")

  def makeId(self):
    """
    Creates an identifier for the reaction to uniquely
    identifies the reactants and products.
    :return str:
    """
    def joinMoleculeNames(molecules):
      names = [m.molecule.name for m in molecules]
      names.sort()
      return ' + '.join(names)
    #
    identifier = "%s %s %s" % (
        joinMoleculeNames(self.reactants),
        REACTION_SEPARATOR,
        joinMoleculeNames(self.products),
        )
    return identifier

  def isEqual(self, other_reaction):
    """
    Checks if two reactions are the same.
    :param Reaction other_reaction:
    :return bool:
    """
    return self.identifier == other_reaction.identifier

  @classmethod
  def initialize(cls, simple):
    """
    :param SimpleSBML simple:
    """
    cls.reactions = []
    [Reaction(r) for r in simple.reactions]
  
