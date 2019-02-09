"""Commonly used utilities."""

from SBMLLint.common import constants as cn

import tellurium as te
import tesbml

TYPE_ANTIMONY = "type_antimony"
TYPE_XML = "type_xml"
TYPE_FILENAME = "type_filename"
XML_HEADER = '<?xml version="1.0" encoding="UTF-8"?>'

def getSBMLDocument(model_reference):
  """
  :param str model_reference:
     sbml xml string
     filename
     antimony source
  :return tesbml.libsbml.Model:
  Returning an SBML document seems to avoid a segementation
  fault from the interations of tesbml and tellurium.
  """
  # Determine the type of the model_reference
  if "<sbml" in model_reference:
    reference_type = TYPE_XML
  elif ("->" in model_reference) or (";" in model_reference):
    reference_type = TYPE_ANTIMONY
  else:
    reference_type = TYPE_FILENAME
  # Process each type
  reader = tesbml.SBMLReader()
  if reference_type == TYPE_FILENAME:
    document = reader.readSBML(model_reference)
  elif reference_type == TYPE_XML:
    document = reader.readSBMLFromString(model_reference)
  else:
    sbml = getSBMLStringFromAntimony(model_reference)
    document = reader.readSBMLFromString(sbml)
  # Extract the model
  if (document.getNumErrors() > 0):
    raise IOError("Errors in SBML document\n%s" 
        % document.printErrors())
  return document

def getSBMLStringFromAntimony(antimony_stg):
  """
  Constructs an SBML model from the antimony string.
  :param str antimony_stg:
  :return str: SBML model in xml format
  """
  rr = te.loada(antimony_stg)
  sbml = rr.getSBML()
  return sbml

def makeSBMLFile(stg, path=cn.TEST_FILE2):
  sbmlstr = getSBMLStringFromAntimony(stg)
  xmlstr = "%s\n%s" % (XML_HEADER, sbmlstr)
  with open(path, 'w') as fd:
    fd.writelines(xmlstr)
