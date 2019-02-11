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
import re
import tesbml
import unittest


IGNORE_TEST = False
UNIUNI = 'SHMTr'
MULTIMULTI = 'GARFT'
NAMEFILTER = "[0-9a-zA-Z]+"


#############################
# Tests
#############################
class TestSOM(unittest.TestCase):

	def setUp(self):
		if IGNORE_TEST:
			return
		self.simple = SimpleSBML(cn.TEST_FILE3)
		Reaction(self.simple._model.getReaction(UNIUNI))
		Reaction(self.simple._model.getReaction(MULTIMULTI))
		self.reaction1 = Reaction.reactions[0]
		self.reaction2 = Reaction.reactions[1]
		self.molecules = Molecule.molecules
		SOM.initialize(self.molecules)
		self.soms = SOM.soms

	def testInitialize(self):
		if IGNORE_TEST:
			return
		SOM.initialize(self.molecules)
		self.assertEqual(len(SOM.soms), len(Molecule.molecules))

	def testMakeId(self):
		if IGNORE_TEST:
			return
		self.assertEqual(len(self.soms), len(self.molecules))
		som = self.soms[0]
		self.assertTrue(som.molecules.intersection(Molecule.molecules))
		molecule = list(re.findall(NAMEFILTER, som.identifier))[0]
		self.assertEqual(list(som.molecules)[0], Molecule.getMolecule(molecule))

	def testMerge(self):
		if IGNORE_TEST:
			return
		SOM.merge(self.reaction1)
		self.assertGreater(len(Molecule.molecules), len(SOM.soms))


	def testFindSOM(self):
		if IGNORE_TEST:
			return
		new_som = SOM.merge(self.reaction1)
		new_som_molecules = list(new_som.molecules)
		self.assertEqual(SOM.findSOM(new_som_molecules[0]), SOM.findSOM(new_som_molecules[1]))


if __name__ == '__main__':
  unittest.main()
