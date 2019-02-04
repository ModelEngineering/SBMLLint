"""Chemical Reaction."""

from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule


class Reaction(object):
  reactions = []

  def __init__(self, libsbml_reaction):
    self._libsbml_reaction = libsbml_reaction
    self.reactants =  self._getMolecules(
        libsbml_reaction.getReactant, 
        libsbml_reaction.getNumReactants)
    self.products=  self._getMolecules(
        libsbml_reaction.getProduct, 
        libsbml_reaction.getNumProducts)
    self.category = self._getCategory()

  def _getMolecules(self, func_getOne, func_getNum):
    """
    Constructs molecules for the species returned by function.
    :param Function func_getOne: a function with an integer argument
        that returns a species
    :param Function func_getNum: a function with no argument that
        returns the total number of species
    :return list-Molecule:
    """
    species = [reaction.func_getOne(n) 
        for n in range(reaction.func_getNum())]
    return [Molecule(s.getId(), species=s) for s in species]

  def _getCategory(self):
    """
    :return str: reaction category
    """
    num_reactants = len(self.reactants)
    num_products = len(self.products)
    category = [r.category for r in cn.REACTION_CATEGORIES
        if r.predicate(num_reactants, num_products)[0]
    return category

  @clasmethod
  def initialize(cls, simple):
    """
    :param SimpleSBML simple:
    """
    [cls.reactions.append(Reaction(r)) for r in simple.reactions]
  
