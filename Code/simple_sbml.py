"""
Provides simplified, read-only access to an SBML model.
1. Create objects
2. Handle stoichiometry
"""
import sys
import os.path
import libsbml


class SimpleSBML(object):
  """
  Provides access to reactions, species, and parameters.
  """

  def __init__(self, filename):
    """
    :param str filename: File containing the SBML document
    :raises IOError: Error encountered reading the SBML document
    """
    self._filename = filename
    reader = libsbml.SBMLReader()
    document = reader.readSBML(self._filename)
    if (document.getNumErrors() > 0):
      raise IOError("Errors in SBML document\n%s" 
          % document.printErrors())
    self._model = document.getModel()
    self._reactions = self._getReactions()
    self._parameters = self._getParameters()

  def _getReactions(self):
    """
    :param libsbml.Model:
    :return list-of-reactions
    """
    num = self._model.getNumReactions()
    return [self._model.getReaction(n) for n in range(num)]

  def _getParameters(self):
    """
    :param libsbml.Model:
    :return list-of-reactions
    """
    return [self._model.getParameter(n) 
         for n in range(self._model.getNumParameters())]

  def getReactions(self):
    return self._reactions

  def getParameters(self):
    return self._parameters

  def getReactants(self, reaction):
    """
    :param libsbml.Reaction:
    :return list-of-libsbml.SpeciesReference:
    """
    return [reaction.getReactant(n) for n in range(reaction.getNumReactants())]

  def getProducts(self, reaction):
    """
    :param libsbml.Reaction:
    :return list-of-libsbml.SpeciesReference:
    """
    return [reaction.getProduct(n) for n in range(reaction.getNumProducts())]

  @classmethod
  def getReactionString(cls, reaction):
    """
    Provides a string representation of the reaction
    :param libsbml.Reaction reaction:
    """
    reaction_str = ''
    base_length = len(reaction_str)
    for reference in getReactants(reaction):
      if len(reaction_str) > base_length:
        reaction_str += " + " + reference.species
      else:
        reaction_str += reference.species
    reaction_str += "-> "
    base_length = len(reaction_str)
    for reference in getProducts(reaction):
      if len(reaction_str) > base_length:
        reaction_str += " + " + reference.species
      else:
        reaction_str += reference.species
    kinetics_terms = cls.getReactionKineticsTerms(reaction)
    reaction_str += "; " + ", ".join(kinetics_terms)
    return reaction_str

  @classmethod
  def getReactionKineticsTerms(cls, reaction):
    """
    Gets the terms used in the kinetics law for the reaction
    :param libsbml.Reaction
    :return list-of-str: names of the terms
    """
    terms = []
    law = reaction.getKineticLaw()
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


if __name__ == '__main__':
  main(sys.argv)  
