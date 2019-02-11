"""
Tests for Set Of Moletules (SOM)
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common.reaction import Reaction
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.games.som import SOM
from SBMLLint.common import simple_sbml

import numpy as np
import os
import tesbml
import unittest


IGNORE_TEST = False
UNIUNI = 'SHMTr'
MULTIMULTI = 'GARFT'


#############################
# Tests
#############################
class TestSOM(unittest.TestCase):

	def setUp(self):
		self.simple = SimpleSBML(cn.TEST_FILE3)
		Reaction(simple._model.getReaction(UNIUNI))
		Reaction(simple._model.getReaction(MULTIMULTI))
		self.reaction1 = Reaction.reactions[0]
		self.reaction2 = Reaction.reactions[1]
		self.molecules = Molecule.molecules
		SOM.soms = []

	def testMakeId(self):
 		SOM.initialize(self.molecules)
 		self.assertEqual(len(SOM.soms), len(self.molecules))


if __name__ == '__main__':
  unittest.main()
