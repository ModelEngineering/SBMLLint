"""Commonly used utilities."""

from SBMLLint.common import constants as cn
from SBMLLint.common.tellurium_sandbox import TelluriumSandbox

import os

TYPE_ANTIMONY = "type_antimony"
TYPE_XML = "type_xml"
TYPE_FILENAME = "type_filename"
XML_HEADER = '<?xml version="1.0" encoding="UTF-8"?>'

def getXML(model_reference):
  """
  :param str model_reference: 
      the input may be a file reference or a model string
          and the file may be an xml file or an antimony file.
      if it is a model string, it may be an xml string or antimony.
  :raises IOError: Error encountered reading the SBML document
  :return str SBML xml"
  """
  # Check for a file path
  if os.path.isfile(model_reference):
    with open(model_reference, 'r') as fd:
      lines = fd.readlines()
    model_str = ''.join(lines)
  else:
    # Must be a string representation of a model
    model_str = model_reference
  # Process model_str into a model  
  if not "<sbml" in model_str:
    # Antimony
    model_str = getXMLFromAntimony(model_str)
  return model_str

def getXMLFromAntimony(antimony_stg):
  """
  Constructs an SBML model from the antimony string.
  :param str antimony_stg:
  :return str: SBML model in xml format
  """
  sandbox = TelluriumSandbox()
  sandbox.run("getSBMLFromAntimony", antimony_stg)
  if sandbox.return_code != 0:
    raise ValueError("Bad antimony string: %s" % antimony_stg)
  return sandbox.output

def isInt(obj):
  try:
    return str(int(obj)) == str(obj)
  except:
    return False

def isFloat(obj):
  try:
    value = float(obj)
  except:
    return False
  return True

def isSBMLModel(obj):
  """
  Tests if object is a libsbml model
  """
  cls_stg = str(type(obj))
  if ('Model' in cls_stg) and ('lib' in cls_stg):
    return True
  else:
    return False

def uniqueify(collection):
  """
  Prunes the collection so that only unique objects are present.
  Elements of the collection must have the method "isEqual" that
  takes as an argument another member of the collection.
  :param list-obj collection
  :return list-obj:
  """
  result = []
  for ele in collection:
    if all([not ele.isEqual(r) for r in result]):
      result.append(ele)
  return result
     
def checkSBMLDocument(document, model_reference=""): 
  if (document.getNumErrors() > 0):
    raise ValueError("Errors in SBML document\n%s" 
        % model_reference)
