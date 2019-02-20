"""
Tests for Mass Equality Set Graph (MESGraph)
"""
from SBMLLint.common import constants as cn
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
    self.simple = SimpleSBML()
    self.simple.initialize(cn.TEST_FILE5)
    SOM.soms = []
    SOM.initialize(self.simple.molecules)
    self.molecules = self.simple.molecules
    self.mesgraph = MESGraph(soms=SOM.soms)
  
  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(len(self.mesgraph.nodes), INITIAL_NODES)
    self.assertEqual(len(self.mesgraph.edges), INITIAL_EDGES)
    dfg = SOM.findSOM(self.simple.getMolecule(DFG))
    # FIXME: Is there supposed to be a test here?
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
    self.mesgraph.processUniUniReaction(
        self.simple.reactions[UNIUNI0])
    dfg = SOM.findSOM(self.simple.getMolecule(DFG))
    e1 = SOM.findSOM(self.simple.getMolecule(E1))
    self.assertTrue(self.mesgraph.has_node(dfg))
    self.assertTrue(self.mesgraph.has_node(e1))
    self.assertEqual(dfg, e1)
  
  def testProcessUniMultiReaction(self):
    if IGNORE_TEST:
      return
    unimulti_reaction = self.simple.reactions[UNIMULTI]
    self.mesgraph.processUniMultiReaction(unimulti_reaction)
    prods = [SOM.findSOM(product.molecule) for product in unimulti_reaction.products]
    dfg = SOM.findSOM(self.simple.getMolecule(DFG))
    for prod in prods:
      self.assertTrue(self.mesgraph.has_edge(prod, dfg))
  
  def testProcessMultiUniReaction(self):
    if IGNORE_TEST:
      return
    multiuni_reaction = self.simple.reactions[MULTIUNI]
    self.mesgraph.processMultiUniReaction(multiuni_reaction)
    reacts = [SOM.findSOM(reactant.molecule) for reactant in multiuni_reaction.reactants]
    mel = SOM.findSOM(self.simple.getMolecule(MEL))
    for react in reacts:
      self.assertTrue(self.mesgraph.has_edge(react, mel))
  
  def testAddArc(self):
    if IGNORE_TEST:
      return
    source = [self.simple.getMolecule(FRU), 
        self.simple.getMolecule(GLY)]
    destination = [self.simple.getMolecule(E2)]
    dummy_reaction = self.simple.reactions[INEQUAL2]
    self.mesgraph.addArc(source, destination, dummy_reaction)
    self.assertTrue(self.mesgraph.has_edge(SOM.findSOM(source[0]), SOM.findSOM(destination[0])))
    self.assertTrue(self.mesgraph.has_edge(SOM.findSOM(source[1]), SOM.findSOM(destination[0])))
  
  def testCheckTypeOneError(self):
    if IGNORE_TEST:
      return
    uniuni_reaction1 = self.simple.reactions[UNIUNI1]
    uniuni_reaction2 = self.simple.reactions[UNIUNI2]
    inequality_reaction1 = self.simple.reactions[INEQUAL1]
    inequality_reaction2 = self.simple.reactions[INEQUAL2]
    self.mesgraph.processUniUniReaction(uniuni_reaction1)
    self.mesgraph.processUniUniReaction(uniuni_reaction2)
    aa = self.simple.getMolecule(AA)
    cn = self.simple.getMolecule(CN)
    mg = self.simple.getMolecule(MG)
    self.assertTrue(self.mesgraph.checkTypeOneError((aa, cn), inequality_reaction1))
    self.assertFalse(self.mesgraph.checkTypeOneError((mg, aa), inequality_reaction2))
  
  def testAnalyze(self):
    if IGNORE_TEST:
      return
    self.mesgraph = MESGraph(SOM.soms)
    self.mesgraph.analyze(self.simple.reactions)
    self.assertEqual(len(self.mesgraph.nodes), FINAL_NODES)
    self.assertEqual(len(self.mesgraph.edges), FINAL_EDGES)

if __name__ == '__main__':
  unittest.main()

