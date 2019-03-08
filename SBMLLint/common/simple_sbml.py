"""
Pythonic representation of an SBML model with some extensions.
- moietys (functional groups within a molecule)
- molecules (species)
- reactions
SimpleSBML extracts all information required from an SBML model to avoid
saving the libsbml object (since these objects are fragile with python
garbage collection).
"""

from SBMLLint.common import constants as cn
from SBMLLint.common.moiety import Moiety, MoietyStoichiometry
from SBMLLint.common.molecule import Molecule, MoleculeStoichiometry
from SBMLLint.common.reaction import Reaction
from SBMLLint.common import util

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

# filename: name of file processed
# number: index of item
# model: libsbml.Model
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

  def initialize(self, model_reference):
    """
    Initializes the instance variables in the model.
    :param str or libsbml.model: for str may be path or model string
       and file/str may be xml or antimony.
    """
    if util.isSBMLModel(model_reference):
      model = model_reference
    else:
      xml = util.getXML(model_reference)
      reader = tesbml.SBMLReader()
      document = reader.readSBMLFromString(xml)
      util.checkSBMLDocument(document)
      model = document.getModel()
    # Do the initializations
    self.reactions = self._getReactions(model)
    self.molecules = self._getMolecules()
    self.moietys = self._getMoietys()

  def _getReactions(self, model):
    reactions = []
    for nn in range(model.getNumReactions()):
      simple_reaction = Reaction(model.getReaction(nn))
      reactions.append(simple_reaction)
    return reactions

  def getReaction(self, label):
    """
    :param str label: label for the reaction
    :return Reaction/None:
    """
    reactions = [r for r in self.reactions if r.label == label]
    if len(reactions) > 1:
      raise ValueError("Two reactions with the same label: %s" %
          label)
    if len(reactions) == 0:
      return None
    return reactions[0]

  def _getMoietys(self):
    """
    Sees if there is a valid moiety structure.
    If not, the molecule is a single moiety.
    """
    moietys = []
    for molecule in self.molecules:
      try:
        new_moietys = ([m_s.moiety 
          for m_s in molecule.getMoietyStoichiometrys()])
      except ValueError:
        new_moietys = [Moiety(molecule.name)]
      moietys.extend(new_moietys)
    return util.uniqueify(moietys)

  def _getMolecules(self):
    """
    :return dict: key is species name, value is species object
    """
    molecules = []
    for reaction in self.reactions:
      molecules.extend(MoleculeStoichiometry.getMolecules(
          reaction.reactants))
      molecules.extend(MoleculeStoichiometry.getMolecules(
          reaction.products))
    return util.uniqueify(molecules)

  def getMolecule(self, name):
    """
    Finds and returns molecule with given name
    Return None if there is no such molecules
    :param str name:
    """
    molecules = [m for m in self.molecules if m.name == name]
    if len(molecules) > 1:
      raise ValueError("Duplicate names in simple.molecules.")
    elif len(molecules) == 1:
      return molecules[0]
    else:
      return None

  def add(self, element):
    """
    Adds an element of the type to its list
    """
    type_list = {
        Moiety: self.moietys,
        Molecule: self.molecules,
        Reaction: self.reactions,
        }
    this_list = type_list[element.__class__]
    appended_list = list(this_list)
    appended_list.append(element)
    new_list = util.uniqueify(appended_list)
    if len(new_list) > len(this_list):
      this_list.append(element)
    

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
