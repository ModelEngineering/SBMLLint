"""
Provides simplified, read-only access to an SBML model.
"""
import sys
import os.path
import tellurium as te
import tesbml

INITIAL_PATH ="http://www.ebi.ac.uk/biomodels-main/download?mid=BIOMD"


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
    self._reader = tesbml.SBMLReader()
    self._document = self._reader.readSBML(self._filename)
    if (self._document.getNumErrors() > 0):
      raise IOError("Errors in SBML document\n%s" 
          % self._document.printErrors())
    self._model = self._document.getModel()
    self._reactions = self._getReactions()
    self._parameters = self._getParameters()  # dict with key=name
    self._species = self._getSpecies()  # dict with key=name

  def _getSpecies(self):
    speciess = {}
    for idx in range(self._model.getNumSpecies()):
      species = self._model.getSpecies(idx)
      speciess[species.getId()] = species
    return speciess

  def _getReactions(self):
    """
    :param tesbml.Model:
    :return list-of-reactions
    """
    num = self._model.getNumReactions()
    return [self._model.getReaction(n) for n in range(num)]

  def _getParameters(self):
    """
    :param tesbml.Model:
    :return list-of-reactions
    """
    parameters = {}
    for idx in range(self._model.getNumParameters()):
      parameter = self._model.getParameter(idx)
      parameters[parameter.getId()] = parameter
    return parameters

  def getReactions(self):
    return self._reactions

  def getParameters(self):
    return self._parameters.keys()

  def getReactants(self, reaction):
    """
    :param tesbml.Reaction:
    :return list-of-tesbml.SpeciesReference:
    To get the species name: SpeciesReference.species
    To get stoichiometry: SpeciesReference.getStoichiometry
    """
    return [reaction.getReactant(n) for n in range(reaction.getNumReactants())]

  def getProducts(self, reaction):
    """
    :param tesbml.Reaction:
    :return list-of-tesbml.SpeciesReference:
    """
    return [reaction.getProduct(n) for n in range(reaction.getNumProducts())]

  @classmethod
  def getReactionString(cls, reaction):
    """
    Provides a string representation of the reaction
    :param tesbml.Reaction reaction:
    """
    reaction_str = ''
    base_length = len(reaction_str)
    for reference in etReactants(reaction):
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
    :param tesbml.Reaction
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

  def isSpecies(self, name):
    """
    Determines if the name is a chemical species
    """
    return self._species.has_key(name)

  def isParameter(self, name):
    """
    Determines if the name is a parameter
    """
    return self._parameters.has_key(name)




def biomodelIterator(initial=1, final=1000):
  """
  :return int, libsbml.model: BioModels number, Model
  """
  num = initial - 1
  for _ in range(final-initial+1):
    num += 1
    formatted_num = format(num, "010")
    url = "%s%s" % (INITIAL_PATH, formatted_num)
    try:
      rr = te.tellurium.RoadRunner(url)
    except:
      break
    model_stg = rr.getCurrentSBML()
    reader = tesbml.SBMLReader()
    document = reader.readSBMLFromString(model_stg)
    yield num, document.getModel()
