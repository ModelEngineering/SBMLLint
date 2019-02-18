"""Chemical Reaction."""

from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule, MoleculeStoichiometry
from SBMLLint.common.simple_sbml import SimpleSBML

import numpy as np


REACTION_SEPARATOR = "->"


################# FUNCTIONS ###################
def getMolecules(libsbml_reaction, func):
  """
  Constructs molecules for the species returned by function.
  :param Function func: gets libsbml.SpeciesReference
  :return list-MoleculeStoichiometry:
  """
  species = func(libsbml_reaction)
  molecule_stoichiometrys = []
  for spc in species:
    molecule = Molecule.getMolecule(spc.species)
    if molecule is None:
        molecule = Molecule(spc.species)
    molecule_stoichiometrys.append(MoleculeStoichiometry(
        molecule,
        spc.getStoichiometry())
        )
  return molecule_stoichiometrys


################# Classes ###################
class Reaction(object):
  reactions = []  # All reactions

  def __init__(self, libsbml_reaction):
    self.reactants = getMolecules(libsbml_reaction,
        SimpleSBML.getReactants)
    self.products = getMolecules(libsbml_reaction,
        SimpleSBML.getProducts)
    self.kinetic_law = libsbml_reaction.getKineticLaw().formula
    self.label = libsbml_reaction.getId()
    self.identifier = self.makeIdentifier(is_include_kinetics=True)
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

  def makeIdentifier(self, is_include_kinetics=True):
    """
    Provides a string representation of the reaction
    :param bool is_include_kinetics: include the kinetics formula
    :return str:
    """
    def makeStoichiometryString(molecule_stoichiometry):
      num = molecule_stoichiometry.stoichiometry
      if np.isclose(num, 1.0):
        return ''
      else:
        return "%2.2f " % num
    #
    def makeTermCollection(molecule_stoichiometries):
      """
      Formats a set of terms with stoichiometries.
      :param list-MoleculeStoichiometry:
      :return str:
      """
      term_collection = ''
      for m_s in molecule_stoichiometries:
        term = "%s%s" % (makeStoichiometryString(m_s), str(m_s.molecule))
        if len(term_collection) == 0:
          term_collection += term
        else:
          term_collection += " + " + term
      return term_collection
    #
    reactant_collection = makeTermCollection(self.reactants)
    product_collection = makeTermCollection(self.products)
    if is_include_kinetics:
      if self.kinetic_law is None:
        formula_str = ''
      else:
        formula_str = "; " + self.kinetic_law
    else:
      formula_str = ''
    reaction_str = "%s: %s -> %s" % (self.label,
        reactant_collection, product_collection)
    reaction_str = reaction_str + formula_str
    return reaction_str

  def __repr__(self):
    return self.getId()

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
  
