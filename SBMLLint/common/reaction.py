"""Chemical Reaction."""

from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule, MoleculeStoichiometry

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

  def __init__(self, libsbml_reaction):
    self.reactants = self.makeMoleculeStoichiometrys(
        libsbml_reaction.getReactant,
        libsbml_reaction.getNumReactants)
    self.products = self.makeMoleculeStoichiometrys(
        libsbml_reaction.getProduct,
        libsbml_reaction.getNumProducts)
    if libsbml_reaction.getKineticLaw() is not None:
      self.kinetics_law = libsbml_reaction.getKineticLaw().formula
    else:
      self.kinetics_law = None
    self.label = libsbml_reaction.getId()
    self.identifier = self.makeIdentifier(is_include_kinetics=True)
    self.category = self.getCategory()
    self.kinetics_terms = self.getKineticsTerms(libsbml_reaction)

  def makeMoleculeStoichiometrys(self, func_get_one, func_get_num):
    """
    Creates a list of MoleculeStoichiometry
    :param Function func_get_one: get one element by index
    :param Function func_get_num: get number of elements
    :return list-MoleculeStoichiometry:
    """
    result = []
    collection = [func_get_one(n) for n in range(func_get_num())]
    for s_r in collection:
      molecule = Molecule(s_r.species)
      stoich = s_r.getStoichiometry()
      result.append(MoleculeStoichiometry(molecule, stoich))
    return result

  def getId(self, is_include_kinetics=True, is_include_label=True):
    """
    Constructs an ID that may be a substring
    of the the full reaction identifier.
    :param bool is_include_kinetics: Include the kinetics law
    :return str:
    """
    result = self.identifier
    if not is_include_kinetics:
      pos = result.find(cn.KINETICS_SEPARATOR)
      if pos > 0:
        result = result[:pos]
    if not is_include_label:
      pos = result.find(cn.LABEL_SEPARATOR)
      if pos > 0:
        result = result[pos+2:]  # Eliminate the separator and space
    return result

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
      if self.kinetics_law is None:
        formula_str = ''
      else:
        formula_str = "; " + self.kinetics_law
    else:
      formula_str = ''
    reaction_str = "%s: %s -> %s" % (self.label,
        reactant_collection, product_collection)
    reaction_str = reaction_str + formula_str
    return reaction_str

  def __repr__(self):
    return self.identifier

  def getCategory(self):
    """
    :return str: reaction category
    """
    num_reactants = len([r.molecule for r in self.reactants \
                         if r.molecule.name!=cn.EMPTYSET])
    num_products = len([p.molecule for p in self.products \
                        if p.molecule.name!=cn.EMPTYSET])
    stoichiometry_reactants = sum([r.stoichiometry for r \
                                   in self.reactants \
                                   if r.molecule.name!=cn.EMPTYSET])
    stoichiometry_products = sum([p.stoichiometry for p \
                                  in self.products \
                                  if p.molecule.name!=cn.EMPTYSET])
    for reaction_category in cn.REACTION_CATEGORIES:
      if reaction_category.predicate(num_reactants, num_products, 
                                     stoichiometry_reactants, 
                                     stoichiometry_products):
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

  def getKineticsTerms(self, libsbml_reaction):
    """
    Gets the terms used in the kinetics law for the reaction
    :param libsbml.Reaction libsbml_reaction:
    :return list-of-str: names of the terms
    """
    terms = []
    law = libsbml_reaction.getKineticLaw()
    if law is not None:
      math = law.getMath()
      asts = [math]
      while len(asts) > 0:
        this_ast = asts.pop()
        if this_ast.isName():
          terms.append(this_ast.getName())
        else:
          pass
        num = this_ast.getNumChildren()
        for idx in range(num):
          asts.append(this_ast.getChild(idx))
    return terms

  @classmethod
  def find(cls, reactions, category=cn.REACTION_1_1):
    return [r for r in reactions if r.category == category]
  
