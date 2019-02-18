"""
Tests for Mass Equality Set Graph (MESGraph)
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common.reaction import Reaction
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.games.mesgraph import MESGraph
from SBMLLint.games.som import SOM
from SBMLLint.common import simple_sbml

import numpy as np
import os
import re
import tesbml
import unittest


IGNORE_TEST = False
INITIAL_NODES = 14
INITIAL_EDGES = 0
FINAL_NODES = 6
FINAL_EDGES = 6
UNIUNI0 = 0
UNIUNI1 = 7
UNIUNI2 = 8
UNIMULTI = 2
MULTIUNI = 13
INEQUAL1 = 14
INEQUAL2 = 15
AA = "AA"
CN = "Cn"
DFG = "DFG"
E1 = "E1"
E2 = "E2"
FRU = "Fru"
GLY = "Gly"
MEL = "Mel"
MG = "MG"

#############################
# Tests
#############################
class TestMESGraph(unittest.TestCase):

	def setUp(self):
		if IGNORE_TEST:
			return
		self.simple = SimpleSBML(cn.TEST_FILE5)
		Reaction.initialize(self.simple)
		SOM.initialize(Molecule.molecules)
		self.molecules = Molecule.molecules
		self.mesgraph = MESGraph(SOM.soms)

	def testConstructor(self):
		if IGNORE_TEST:
			return
		self.assertEqual(len(self.mesgraph.nodes), INITIAL_NODES)
		self.assertEqual(len(self.mesgraph.edges), INITIAL_EDGES)
		dfg = SOM.findSOM(Molecule.getMolecule(DFG))
		self.mesgraph.has_node(dfg)

	def testMakeId(self):	
		if IGNORE_TEST:
			return
		identifier = ""
		for key, som in enumerate(SOM.soms):
		    identifier = identifier + som.identifier
		    if key < len(SOM.soms)-1:
		        identifier = identifier + ";"
		self.assertEqual(identifier, self.mesgraph.identifier)

	def testProcessUniUniReaction(self):
		if IGNORE_TEST:
			return
		self.mesgraph.processUniUniReaction(Reaction.reactions[UNIUNI0])
		dfg = SOM.findSOM(Molecule.getMolecule(DFG))
		e1 = SOM.findSOM(Molecule.getMolecule(E1))
		self.assertTrue(self.mesgraph.has_node(dfg))
		self.assertTrue(self.mesgraph.has_node(e1))
		self.assertEqual(dfg, e1)

	def testProcessUniMultiReaction(self):
		if IGNORE_TEST:
			return
		unimulti_reaction = Reaction.reactions[UNIMULTI]
		self.mesgraph.processUniMultiReaction(unimulti_reaction)
		prods = [SOM.findSOM(product.molecule) for product in unimulti_reaction.products]
		dfg = SOM.findSOM(Molecule.getMolecule(DFG))
		for prod in prods:
		    self.assertTrue(self.mesgraph.has_edge(prod, dfg))

	def testProcessMultiUniReaction(self):
		if IGNORE_TEST:
			return
		multiuni_reaction = Reaction.reactions[MULTIUNI]
		self.mesgraph.processMultiUniReaction(multiuni_reaction)
		reacts = [SOM.findSOM(reactant.molecule) for reactant in multiuni_reaction.reactants]
		mel = SOM.findSOM(Molecule.getMolecule(MEL))
		for react in reacts:
		    self.assertTrue(self.mesgraph.has_edge(react, mel))

	def testAddArc(self):
		if IGNORE_TEST:
			return
		source = [Molecule.getMolecule(FRU), Molecule.getMolecule(GLY)]
		destination = [Molecule.getMolecule(E2)]
		dummy_reaction = Reaction.reactions[INEQUAL2]
		self.mesgraph.addArc(source, destination, dummy_reaction)
		self.assertTrue(self.mesgraph.has_edge(SOM.findSOM(source[0]), SOM.findSOM(destination[0])))
		self.assertTrue(self.mesgraph.has_edge(SOM.findSOM(source[1]), SOM.findSOM(destination[0])))

	def testCheckTypeOneError(self):
		if IGNORE_TEST:
			return
		uniuni_reaction1 = Reaction.reactions[UNIUNI1]
		uniuni_reaction2 = Reaction.reactions[UNIUNI2]
		inequality_reaction1 = Reaction.reactions[INEQUAL1]
		inequality_reaction2 = Reaction.reactions[INEQUAL2]
		self.mesgraph.processUniUniReaction(uniuni_reaction1)
		self.mesgraph.processUniUniReaction(uniuni_reaction2)
		aa = Molecule.getMolecule(AA)
		cn = Molecule.getMolecule(CN)
		mg = Molecule.getMolecule(MG)
		self.assertTrue(self.mesgraph.checkTypeOneError((aa, cn), inequality_reaction1))
		self.assertFalse(self.mesgraph.checkTypeOneError((mg, aa), inequality_reaction2))

	def testAnalyze(self):
		if IGNORE_TEST:
			return
		self.mesgraph = MESGraph(SOM.soms)
		self.mesgraph.analyze(Reaction.reactions)
		self.assertEqual(len(self.mesgraph.nodes), FINAL_NODES)
		self.assertEqual(len(self.mesgraph.edges), FINAL_EDGES)
		
if __name__ == '__main__':
  unittest.main()

