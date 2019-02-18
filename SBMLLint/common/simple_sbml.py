"""
Pythonic representation of an SBML model with some extensions.
- moietys (functional groups within a molecule)
- molecules (species)
- reactions
"""

from SBMLLint.common import constants as cn
import collections
import os.path
import numpy as np
import sys
import tesbml
import urllib3
import warnings


TYPE_MODEL = "type_model"  # libsbml model
TYPE_XML = "type_xml"  # XML string
TYPE_ANTIMONY = "type_xml"  # Antimony string
TYPE_FILE = "type_file" # File reference

IteratorItem = collections.namedtuple('IteratorItem',
    'filename number model')


class SimpleSBML(object):
  """
  This class address stability of the underlying tesbml 
  (libsbml) library that seems not to survive garbage collection
  by python (e.g., returning a libsbml object to a caller.) As
  a result, no libsbml object is maintained by SimpleSBML instances.
  """

  def __init__(self):
    """
    Initializes instance variables
    """
    self.moietys = []
    self.molecules = []
    self.reactions = []

  def _getSBMLModel(self, model_reference):
    """
    :param str or tesbml.libsbml.model model_reference: 
        if the input is a str it may be a file reference or a model string
            and the file may be an xml file or an antimony file.
        if it is a model string, it may be an xml string or antimony.
    :raises IOError: Error encountered reading the SBML document
    """
    # check for a libsbml model
    if 'Model' in str(type(model_reference)):
      return model_reference
    if not isinstance(model_reference, str):
      raise ValueError("Invalid model_reference: %s" % str(model_reference))
    # Check for a file path
    if os.path.isfile(model_reference):
      with open(model_reference, 'r') as fd:
        lines = fd.readlines()
      model_str = ''.join(lines)
    else:
      # Must be a string representation of a model
      model_str = model_reference
    # Process model_str into a model  
    if "<sbml" in model_str:
      # SBML file
      reader = tesbml.SBMLReader()
      document = reader.readSBMLFromString(model_str)
      if (document.getNumErrors() > 0):
        raise ValueError("Errors in SBML document\n%s" 
            % model_reference)
      model = self._document.getModel()
    else:
      # Assume this is an antimony string
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

  def getSpeciesNames(self):
    """
    :return list-str: list of species names
    """
    return list(self.species.keys())

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
