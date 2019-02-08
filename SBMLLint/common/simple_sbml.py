"""
Provides simplified, read-only access to an SBML model.
"""

from SBMLLint.common import constants as cn
import collections
import os.path
import numpy as np
import sys
import tesbml
import urllib3
import warnings


TYPE_MODEL = "type_model"
TYPE_XML = "type_xml"
TYPE_FILE = "type_file"

IteratorItem = collections.namedtuple('IteratorItem',
    'filename number model')


class SimpleSBML(object):
  """
  Provides access to reactions, species, and parameters.
  """

  def __init__(self, model_reference):
    """
    :param str or tesbml.libsbml.model model_reference: 
        File  or sbml model
    :raises IOError: Error encountered reading the SBML document
    """
    if isinstance(model_reference, str):
      if "<sbml" in model_reference:
        model_type = TYPE_XML
      else:
        model_type = TYPE_FILE
    else:
      model_type = TYPE_MODEL
    if model_type == TYPE_FILE:
      self._filename = model_reference
      self._reader = tesbml.SBMLReader()
      self._document = self._reader.readSBML(self._filename)
      if (self._document.getNumErrors() > 0):
        raise IOError("Errors in SBML document\n%s" 
            % self._document.printErrors())
      self._model = self._document.getModel()
    elif model_type == TYPE_XML:
      self._filename = None
      self._reader = tesbml.SBMLReader()
      self._document = self._reader.readSBMLFromString(
          model_reference)
      if (self._document.getNumErrors() > 0):
        raise IOError("Errors in SBML document\n%s" 
            % self._document.printErrors())
      self._model = self._document.getModel()
    else:
      self._filename = None
      self._reader = None
      self._document = None
      self._model = model_reference
    self.reactions = self._getReactions()
    self.parameters = self._getParameters()  # dict with key=name
    self.species = self._getSpecies()  # dict with key=name

  def _getReactions(self):
    """
    :param tesbml.libsbml.Model:
    :return list-of-reactions
    """
    num = self._model.getNumReactions()
    return [self._model.getReaction(n) for n in range(num)]

  def _getSpecies(self):
    """
    :return dict: key is species name, value is species object
    """
    speciess = {}
    for idx in range(self._model.getNumSpecies()):
      species = self._model.getSpecies(idx)
      speciess[species.getId()] = species
    return speciess

  def _getReactions(self):
    """
    :param tesbml.libsbml.Model:
    :return list-of-reactions
    """
    num = self._model.getNumReactions()
    return [self._model.getReaction(n) for n in range(num)]

  def _getParameters(self):
    """
    :param tesbml.libsbml.Model:
    :return list-of-reactions
    """
    parameters = {}
    for idx in range(self._model.getNumParameters()):
      parameter = self._model.getParameter(idx)
      parameters[parameter.getId()] = parameter
    return parameters

  def getReactions(self):
    return self.reactions

  def getParameters(self):
    return self.parameters.keys()

  @staticmethod
  def getReactants(reaction):
    """
    :param tesbml.libsbml.Reaction:
    :return list-of-libsbml.SpeciesReference:
    To get the species name: SpeciesReference.species
    To get stoichiometry: SpeciesReference.getStoichiometry
    """
    return [reaction.getReactant(n) for n in range(reaction.getNumReactants())]

  @staticmethod
  def getProducts(reaction):
    """
    :param tesbml.libsbml.Reaction:
    :return list-of-tesbml.libsbml.SpeciesReference:
    """
    return [reaction.getProduct(n) for n in range(reaction.getNumProducts())]

  @classmethod
  def getReactionString(cls, reaction, is_include_kinetics=True):
    """
    Provides a string representation of the reaction
    :param tesbml.libsbml.Reaction reaction:
    :param bool is_include_kinetics: include the kinetics formula
    :return str:
    """
    def makeStoichiometryString(species_reference):
      num = species_reference.getStoichiometry()
      if np.isclose(num, 1.0):
        return ''
      else:
        return "%2.2f " % num
    #
    def makeTermCollection(species_references):
      """
      Formats a set of terms with stoichiometries.
      :param list-SpeciesReference:
      :return str:
      """
      term_collection = ''
      for reference in species_references:
        term = "%s%s" % (makeStoichiometryString(reference),
          reference.species)
        if len(term_collection) == 0:
          term_collection += term
        else:
          term_collection += " + " + term
      return term_collection
    #
    reactant_collection = makeTermCollection(cls.getReactants(reaction))
    product_collection = makeTermCollection(cls.getProducts(reaction))
    if is_include_kinetics:
      formula_str = "; " + reaction.getKineticLaw().formula
    else:
      formula_str = ''
    reaction_str = "%s: %s -> %s%s" % (reaction.id, reactant_collection,
        product_collection, formula_str)
    return reaction_str

  @staticmethod
  def getReactionKineticsTerms(reaction):
    """
    Gets the terms used in the kinetics law for the reaction
    :param tesbml.libsbml.Reaction
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
    return name in list(self.species.keys())

  def isParameter(self, name):
    """
    Determines if the name is a parameter
    """
    return name in self.parameters.keys()

###################### FUNCTIONS #############################
def readURL(url):
  """
  :param str url:
  :return str: file content
  """
  def do():
    http = urllib3.PoolManager()
    response = http.request('GET', url)
    return response.data.decode("utf-8") 
 # Catch bogus warnings
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    result = do()
  return result
  
def modelIterator(initial=0, final=1000, data_dir=cn.DATA_DIR):
  """
  Iterates across all models in a data directory.
  :param int initial: initial file to process
  :param int final: final file to process
  :param str data_dir: absolute path of the directory containing
      the xml files
  :return IteratorItem:
  """
  files = [f for f in os.listdir(data_dir) if f[-4:] == ".xml"]
  begin_num = max(initial, 0)
  num = begin_num - 1
  end_num = min(len(files), final)
  for filename in files[begin_num:end_num]:
    path = os.path.join(data_dir, filename)
    num += 1
    with open(path, 'r') as fd:
      lines = ''.join(fd.readlines())
      reader = tesbml.libsbml.SBMLReader()
      document = reader.readSBMLFromString(lines)
      model = document.getModel()
    iterator_item = IteratorItem(filename=filename,
        model=model, number=num)
    yield iterator_item
