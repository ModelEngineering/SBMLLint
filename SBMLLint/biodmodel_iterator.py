"""Generator to iterate through biomodels."""

import tellurium as te
import tesbml

INITIAL_PATH ="http://www.ebi.ac.uk/biomodels-main/download?mid=BIOMD"


def biomodel_iterator(initial=1, final=1000):
  """
  :return int, libsbml.model: BioModels number, Model
  """
  num = initial
  for _ in range(final-initial):
    formatted_num = format(num, "010")
    url = "%s%s" % (INITIAL_PATH, formatted_num)
    try:
      rr = te.tellurium.RoadRunner(url)
    except:
      break
    model_stg = rr.getCurrentSBML()
    reader = tesbml.SBMLReader()
    document = reader.readSBMLFromString(stg)
    yield document.getModel()
