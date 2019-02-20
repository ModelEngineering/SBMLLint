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
		self.simple = SimpleSBML()
		self.simple.initialize(cn.TEST_FILE3)
                self.simple.add(Reaction(
		    self.simple._model.getReaction(UNIUNI)))
                self.simple.add(Reaction(
		    self.simple._model.getReaction(MULTIMULTI)))
                self.simple.add(Reaction
		    self.simple._model.getReaction(MULTIUNI)))
		self.uniuni = self.simple.reactions[0]
		self.multimulti = self.simple.reactions[1]
		self.multiuni = self.simple.reactions[2]
		self.molecules = self.simple.molecules
		SOM.initialize(self.molecules)
		self.soms = SOM.soms

	def testInitialize(self):
		if IGNORE_TEST:
			return
		SOM.initialize(self.molecules)
		self.assertEqual(len(SOM.soms), len(self.simple.molecules))

	def testMakeId(self):
		if IGNORE_TEST:
			return
		self.assertEqual(len(self.soms), len(self.molecules))
		som = self.soms[0]
		self.assertTrue(som.molecules.intersection(self.simple.molecules))
		molecule = list(re.findall(NAMEFILTER, som.identifier))[0]
		self.assertEqual(list(som.molecules)[0], self.simple.getMolecule(molecule))

	def testMerge(self):
		if IGNORE_TEST:
			return
		SOM.merge(self.uniuni)
		self.assertGreater(len(self.simple.molecules), len(SOM.soms))

	def testFindSOM(self):
		if IGNORE_TEST:
			return
		simple = SimpleSBML()
		simple.initialize(cn.TEST_FILE3)
		SOM.soms = []
                self.simple.add(Reaction(
		    self.simple._model.getReaction(UNIUNI))
                self.simple.add(Reaction(
		    self.simple._model.getReaction(MULTIMULTI))
                self.simple.add(Reaction(
		    self.simple._model.getReaction(MULTIUNI))
		molecule = simple.getMolecule(MOLECULE)
		SOM.initialize(simple.molecules)
		self.assertEqual(list(SOM.findSOM(molecule).molecules)[0], molecule)
		reactions = Reaction.find(simple.reactions,
                     category=cn.REACTION_1_1)
		new_som = SOM.merge(reactions[0])
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
		reduced_reaction = SOM.reduce(self.multimulti)
		self.assertIsInstance(reduced_reaction, Reaction)
		self.assertEqual(reduced_reaction.category, cn.REACTION_n_n)
		self.assertGreater(num_reactants, len(reduced_reaction.reactants))
		self.assertGreater(num_products, len(reduced_reaction.products))

	def testAddSOM(self):
		if IGNORE_TEST:
			return
		num_soms = len(SOM.soms)
		som_popped = SOM.soms.pop()
		som_existing = SOM.soms[0]

		self.assertEqual(num_soms-1, len(SOM.soms))
		self.assertTrue(som_existing in SOM.soms)
		self.assertFalse(som_popped in SOM.soms)
		SOM.addSOM(som_existing)
		self.assertEqual(num_soms-1, len(SOM.soms))
		SOM.addSOM(som_popped)
		self.assertEqual(num_soms, len(SOM.soms))
		self.assertTrue(som_popped in SOM.soms)

if __name__ == '__main__':
  unittest.main()
