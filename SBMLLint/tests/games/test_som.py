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
UNIUNI = 'R12'
MULTIMULTI = 'R1'
MULTIUNI = 'R13'
MOLECULE = 'halfglucose'
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
		Reaction(self.simple._model.getReaction(MULTIUNI))
		self.uniuni = Reaction.reactions[0]
		self.multimulti = Reaction.reactions[1]
		self.multiuni = Reaction.reactions[2]
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
		SOM.merge(self.uniuni)
		self.assertGreater(len(Molecule.molecules), len(SOM.soms))

	def testFindSOM(self):
		if IGNORE_TEST:
			return
		molecule = Molecule.getMolecule(MOLECULE)
		self.assertEqual(list(SOM.findSOM(molecule).molecules)[0], molecule)
		new_som = SOM.merge(self.uniuni)
		new_som_molecules = list(new_som.molecules)
		self.assertEqual(SOM.findSOM(new_som_molecules[0]), SOM.findSOM(new_som_molecules[1]))

	def testReduce(self):
		if IGNORE_TEST:
			return
		SOM.merge(self.uniuni)
		self.assertFalse(SOM.reduce(self.uniuni))
		self.assertFalse(SOM.reduce(self.multiuni))
		num_reactants = len(self.multimulti.reactants)
		num_products = len(self.multimulti.products)
		reduced_reaction = SOM.reduce(self.multimulti)
		self.assertIsInstance(reduced_reaction, Reaction)
		self.assertEqual(reduced_reaction.category, cn.REACTION_n_n)
		self.assertGreater(num_reactants, len(reduced_reaction.reactants))
		self.assertGreater(num_products, len(reduced_reaction.products))

	def testaddSOM(self):
		if IGNORE_TEST:
			return
		num_soms = len(SOM.soms)
		som_popped = SOM.soms.pop()
		som_existing = SOM.soms[0]

		self.assertEqual(num_soms-1, len(SOM.soms))
		SOM.addSOM(som_existing)
		self.assertEqual(num_soms-1, len(SOM.soms))
		SOM.addSOM(som_popped)
		self.assertEqual(num_soms, len(SOM.soms))

if __name__ == '__main__':
  unittest.main()
