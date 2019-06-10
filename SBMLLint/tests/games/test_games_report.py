"""
Test for Games Report
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common.reaction import Reaction
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.games.games_pp import SOMStoichiometry, SOMReaction, GAMES_PP
from SBMLLint.games.games_report import GAMESReport, NULL_STR, NUM_STAR
from SBMLLint.games.som import SOM
from SBMLLint.common import simple_sbml

import numpy as np
import pandas as pd
import os
import re
import tesbml
import unittest

IGNORE_TEST = False


#############################
# Tests
#############################
class TestGAMESReport(unittest.TestCase):

  def setUp(self):
    self.simple1 = SimpleSBML()
    self.simple1.initialize(cn.TEST_FILE_GAMESREPORT1)

  # def testReportCancelingError(self):
  # 	m = GAMES_PP(self.simple1)
  # 	m.analyze(error_details=False)
  # 	gr = GAMESReport(m)
  # 	canceling_error_report = gr.reportCancelingError(m.canceling_errors)
  # 	correct_report = NULL_STR
  # 	correct_report = correct_report + "\nThe following reaction:\n"
  # 	correct_report = correct_report + "OxidativePhosphorylation: 6.00 ADP + CTtis -> 6.00 ATP\n\n"
  # 	correct_report = correct_report + "has a mass balance error\nbecause of the following equality."
  # 	correct_report = correct_report + "\n\nADP = ATP by reaction(s):\n"
  # 	correct_report = correct_report + "ATPase: ATP -> ADP\n\n"
  # 	correct_report = correct_report + "*"*NUM_STAR + "\n" + "-"*NUM_STAR + "\n"
  # 	self.assertEqual(correct_report, canceling_error_report)











