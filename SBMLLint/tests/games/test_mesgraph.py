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
# Constants for simple
AA = "AA"
CN = "Cn"
DFG = "DFG"
E1 = "E1"
E2 = "E2"
FRU = "Fru"
GLY = "Gly"
MEL = "Mel"
MG = "MG"
# Constants for simple2
REACTION1 = "SHMTr"
REACTION2 = "MTHFR"
REACTION3 = "MTR"


#############################
# Tests
#############################
class TestMESGraph(unittest.TestCase):

  def setUp(self):
    # SimpleSBML with type I error
    self.simple = SimpleSBML()
    self.simple.initialize(cn.TEST_FILE6)
    # SimpleSBML with type II error
    self.simple2 = SimpleSBML()
    self.simple2.initialize(cn.TEST_FILE7)
    self.mesgraph = MESGraph(self.simple)
  
  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(len(self.mesgraph.nodes), INITIAL_NODES)
    self.assertEqual(len(self.mesgraph.edges), INITIAL_EDGES)
    dfg = self.simple.getMolecule(DFG)
    # molecules is a list of one-molecule sets
    molecules = [som.molecules for som in self.mesgraph.nodes]
    self.assertTrue({self.simple.getMolecule(DFG)} in molecules)

  def testInitializeSOMs(self):
    if IGNORE_TEST:
      return
    for node in self.mesgraph.nodes:
      self.assertEqual(type(node), SOM)
  
  def testMakeId(self):
    if IGNORE_TEST:
      return
    identifier = ""
    for key, som in enumerate(self.mesgraph.nodes):
      identifier = identifier + som.identifier
      if key < len(self.mesgraph.nodes)-1:
        identifier = identifier + ";"
    self.assertEqual(identifier, self.mesgraph.identifier)

  def testGetNode(self):
    if IGNORE_TEST:
      return
    aa = self.simple.getMolecule(AA)
    aa_node = self.mesgraph.getNode(aa)
    self.assertEqual(type(aa_node), SOM)
    self.assertEqual(aa_node.molecules, {aa})
  
  def testProcessUniUniReaction(self):
    if IGNORE_TEST:
      return
    self.mesgraph.processUniUniReaction(
        self.simple.reactions[UNIUNI0])
    dfg = self.mesgraph.getNode(self.simple.getMolecule(DFG))
    e1 = self.mesgraph.getNode(self.simple.getMolecule(E1))
    self.assertTrue(self.mesgraph.has_node(dfg))
    self.assertTrue(self.mesgraph.has_node(e1))
    self.assertEqual(dfg, e1)
  
  def testProcessUniMultiReaction(self):
    if IGNORE_TEST:
      return
    unimulti_reaction = self.simple.reactions[UNIMULTI]
    self.mesgraph.processUniMultiReaction(unimulti_reaction)
    prods = [self.mesgraph.getNode(product.molecule) for product in unimulti_reaction.products]
    dfg = self.mesgraph.getNode(self.simple.getMolecule(DFG))
    for prod in prods:
      self.assertTrue(self.mesgraph.has_edge(prod, dfg))
  
  def testProcessMultiUniReaction(self):
    if IGNORE_TEST:
      return
    multiuni_reaction = self.simple.reactions[MULTIUNI]
    self.mesgraph.processMultiUniReaction(multiuni_reaction)
    reacts = [self.mesgraph.getNode(reactant.molecule) for reactant in multiuni_reaction.reactants]
    mel = self.mesgraph.getNode(self.simple.getMolecule(MEL))
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
    arc1 = [self.mesgraph.getNode(source[0]), 
            self.mesgraph.getNode(destination[0])]
    arc2 = [self.mesgraph.getNode(source[1]), 
            self.mesgraph.getNode(destination[0])]
    self.assertTrue(self.mesgraph.has_edge(arc1[0], arc1[1]))
    self.assertTrue(self.mesgraph.has_edge(arc2[0], arc2[1]))
    reaction_label1 = self.mesgraph.get_edge_data(arc1[0], arc1[1])[cn.REACTION][0]
    reaction_label2 = self.mesgraph.get_edge_data(arc2[0], arc1[1])[cn.REACTION][0]
    self.assertEqual(reaction_label1, dummy_reaction.label)
    self.assertEqual(reaction_label1, reaction_label2)
  
  def testGetSOMPath(self):
    pass

  def testPrintSOMPath(self):
    

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
    self.assertTrue(len(self.mesgraph.type_one_errors)>0)
    self.assertFalse(len(self.mesgraph.type_two_errors)>0)

  def testCheckTypeTwoError(self):
    # TODO: need to fix according to stoichiometry; 
    # Need to use another model once Reaction class is fixed
    if IGNORE_TEST:
      return
    mesgraph2 = MESGraph(self.simple2)
    mesgraph2.processUniUniReaction(self.simple2.getReaction(REACTION1))
    mesgraph2.processMultiUniReaction(self.simple2.getReaction(REACTION2))
    mesgraph2.processMultiUniReaction(self.simple2.getReaction(REACTION3))
    self.assertTrue(mesgraph2.checkTypeTwoError())
    self.assertFalse(len(mesgraph2.type_one_errors)>0)
    self.assertTrue(len(mesgraph2.type_two_errors)>0)

  def testAnalyze(self):
    if IGNORE_TEST:
      return
    mesgraph1 = MESGraph(self.simple)
    mesgraph1.analyze(self.simple.reactions)
    self.assertEqual(len(mesgraph1.nodes), FINAL_NODES)
    self.assertEqual(len(mesgraph1.edges), FINAL_EDGES)
    self.assertTrue(len(mesgraph1.type_one_errors)>0)
    self.assertFalse(len(mesgraph1.type_two_errors)>0)
    #
    mesgraph2 = MESGraph(self.simple2)
    mesgraph2.analyze(self.simple2.reactions)
    self.assertTrue(len(mesgraph2.type_one_errors)>0)
    self.assertTrue(len(mesgraph2.type_two_errors)>0)

if __name__ == '__main__':
  unittest.main()

