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
import unittest

IGNORE_TEST = False
# Constants for simple and simple2, non multi-multi reactions
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
# Constants for simple3, Type III errors
V1 = "v1"
V2 = "v2"
V3 = "v3"
AMP = "AMP"
ADP = "ADP"
ATP = "ATP"
# Constants for simple4, Tyep IV errors, error from reducing
ATPASE = "ATPase"
LOWER = "lower"
FRU16P2 = "Fru16P2"
#simple5, type V error, SOM cycle
DIH = "DIH"
DIN = "DIN"
#simple6, no errors in multi-multi model
R5 = "R5"
R12 = "R12"
AC = "Ac"
ACP = "AcP"


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
    # SimpleSBML for type III error
    self.simple3 = SimpleSBML()
    self.simple3.initialize(cn.TEST_FILE10)
    # SimpleSBML for type IV error
    self.simple4 = SimpleSBML()
    self.simple4.initialize(cn.TEST_FILE11)
    # simple5 for type V error    
    self.simple5 = SimpleSBML()
    self.simple5.initialize(cn.TEST_FILE12)
    # simple6 will test multi-multi model with no errors
    self.simple6 = SimpleSBML()
    self.simple6.initialize(cn.TEST_FILE3)
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
    self.assertEqual(len(self.mesgraph.type_one_errors), 0)
    self.assertEqual(len(self.mesgraph.type_two_errors), 0)
    self.assertEqual(len(self.mesgraph.type_three_errors), 0)
    self.assertEqual(len(self.mesgraph.type_four_errors), 0)
    self.assertEqual(len(self.mesgraph.type_five_errors), 0)

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
  
  def testMergeNodes(self):
    if IGNORE_TEST:
      return
    m3 = MESGraph(self.simple3)
    v1_reaction = self.simple3.getReaction(V1)
    v2_reaction = self.simple3.getReaction(V2)
    amp_molecule = self.simple3.getMolecule(AMP)
    adp_molecule = self.simple3.getMolecule(ADP)
    atp_molecule = self.simple3.getMolecule(ATP)
    amp_som = m3.getNode(amp_molecule)
    adp_som = m3.getNode(adp_molecule)
    atp_som = m3.getNode(atp_molecule)
    ampatp_som = m3.mergeNodes(amp_som, atp_som, v1_reaction)
    self.assertTrue(amp_molecule in ampatp_som.molecules)
    self.assertTrue(atp_molecule in ampatp_som.molecules)
    self.assertTrue(ampatp_som in m3.nodes)
    # tests if anohter merge will remove previous nodes
    ampadpatp_som = m3.mergeNodes(ampatp_som, adp_som, v2_reaction)
    self.assertFalse(ampatp_som in m3.nodes)
    self.assertFalse(adp_som in m3.nodes)  
    self.assertTrue(ampadpatp_som in m3.nodes)

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
  
  def testAddMultiMultiReaction(self):
    if IGNORE_TEST:
      return
    m3 = MESGraph(self.simple3)
    v2 = self.simple3.getReaction(V2)
    m3.addMultiMultiReaction(v2)
    self.assertTrue(v2 in m3.multimulti_reactions)

  def testAddTypeThreeError(self):
    if IGNORE_TEST:
      return
    m3 = MESGraph(self.simple3)
    v1 = self.simple3.getReaction(V1)
    v2 = self.simple3.getReaction(V2)
    v3 = self.simple3.getReaction(V3)
    m3.processUniUniReaction(v1)
    m3.processUniMultiReaction(v3)
    amp = m3.getNode(self.simple3.getMolecule(AMP))
    atp = m3.getNode(self.simple3.getMolecule(ATP))
    self.assertTrue(m3.addTypeThreeError(amp, atp, v2))
    self.assertTrue(len(m3.type_three_errors), 1)
    error = m3.type_three_errors[0]
    self.assertEqual(error.node1, amp)
    self.assertEqual(error.node2, atp)
    self.assertEqual(error.reactions, [v2.label])

  def testCheckTypeThreeError(self):
    if IGNORE_TEST:
      return
    m3 = MESGraph(self.simple3)
    v1 = self.simple3.getReaction(V1)
    v2 = self.simple3.getReaction(V2)
    v3 = self.simple3.getReaction(V3)
    m3.processUniUniReaction(v1)
    adp = m3.getNode(self.simple3.getMolecule(ADP))
    atp = m3.getNode(self.simple3.getMolecule(ATP))
    self.assertFalse(m3.checkTypeThreeError(adp, atp, v3))
    m3.processUniMultiReaction(v3) 
    self.assertTrue(m3.checkTypeThreeError(adp, atp, v2))
    self.assertTrue(m3.checkTypeThreeError(atp, adp, v2))

  def testReduceReaction(self):
    if IGNORE_TEST:
      return
    m4 = MESGraph(self.simple4)
    atpase = m4.simple.getReaction(ATPASE)
    lower = m4.simple.getReaction(LOWER)
    self.assertFalse(m4.reduceReaction(atpase))
    m4.processUniUniReaction(atpase)
    reduced_reaction = m4.reduceReaction(lower)
    self.assertEqual(type(reduced_reaction), Reaction)
    self.assertEqual(len(reduced_reaction.reactants), 1)
    self.assertEqual(reduced_reaction.products, [])   
    self.assertEqual(reduced_reaction.reactants[0].molecule.name, FRU16P2)

  def testProcessMultiMultiReactions(self):
    if IGNORE_TEST:
      return
    # type III error
    m3 = MESGraph(self.simple3)
    v1 = self.simple3.getReaction(V1)
    v2 = self.simple3.getReaction(V2)
    v3 = self.simple3.getReaction(V3)
    m3.processUniUniReaction(v1)
    m3.processUniMultiReaction(v3) 
    self.assertEqual(m3.type_three_errors, [])
    self.assertTrue(m3.processMultiMultiReaction(v2))
    self.assertEqual(len(m3.type_three_errors), 1) 
    # type IV error
    m4 = MESGraph(self.simple4)
    atpase = m4.simple.getReaction(ATPASE)
    lower = m4.simple.getReaction(LOWER)
    m4.processUniUniReaction(atpase)  
    self.assertEqual(m4.type_four_errors, []) 
    self.assertTrue(m4.processMultiMultiReaction(lower))
    self.assertEqual(len(m4.type_four_errors), 1)
    # no error
    m6 = MESGraph(self.simple6)
    r5 = m6.simple.getReaction(R5)
    r12 = m6.simple.getReaction(R12)
    m6.processUniUniReaction(r12)
    ac = m6.getNode(self.simple6.getMolecule(AC))
    acp = m6.getNode(self.simple6.getMolecule(ACP))
    self.assertFalse(ac==acp)
    self.assertTrue(m6.processMultiMultiReaction(r5))
    ac = m6.getNode(self.simple6.getMolecule(AC))
    acp = m6.getNode(self.simple6.getMolecule(ACP))
    self.assertTrue(ac==acp)

  def testAddArc(self):
    if IGNORE_TEST:
      return
    source = [self.mesgraph.getNode(self.simple.getMolecule(FRU)), 
              self.mesgraph.getNode(self.simple.getMolecule(GLY))]
    destination = [self.mesgraph.getNode(self.simple.getMolecule(E2))]
    dummy_reaction = self.simple.reactions[INEQUAL2]
    self.mesgraph.addArc(source[0], destination[0], dummy_reaction)
    self.mesgraph.addArc(source[1], destination[0], dummy_reaction)
    arc1 = [source[0], destination[0]]
    arc2 = [source[1], destination[0]]
    self.assertTrue(self.mesgraph.has_edge(arc1[0], arc1[1]))
    self.assertTrue(self.mesgraph.has_edge(arc2[0], arc2[1]))
    reaction_label1 = self.mesgraph.get_edge_data(arc1[0], arc1[1])[cn.REACTION][0]
    reaction_label2 = self.mesgraph.get_edge_data(arc2[0], arc1[1])[cn.REACTION][0]
    self.assertEqual(reaction_label1, dummy_reaction.label)
    self.assertEqual(reaction_label1, reaction_label2)
  
  def testGetSOMPath(self):
    if IGNORE_TEST:
      return
    uniuni_reaction = self.simple.reactions[UNIUNI0]
    self.mesgraph.processUniUniReaction(uniuni_reaction)
    dfg = self.simple.getMolecule(DFG)
    e1 = self.simple.getMolecule(E1)
    som = self.mesgraph.getNode(dfg)
    som_path = self.mesgraph.getSOMPath(som, e1, dfg)
    self.assertEqual(type(som_path[0]), cn.PathComponents)
    self.assertEqual(som_path[0].reactions, [uniuni_reaction.label])
    self.assertEqual(len(som_path), 1)

  def testPrintSOMPath(self):
    if IGNORE_TEST:
      return
    uniuni_reaction = self.simple.reactions[UNIUNI0]
    self.mesgraph.processUniUniReaction(uniuni_reaction)
    self.assertTrue(type(self.mesgraph.printSOMPath(DFG, E1))==str)
    self.assertFalse(self.mesgraph.printSOMPath(GLY, MEL))

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
    #
    mesgraph3 = MESGraph(self.simple3)
    mesgraph3.analyze(self.simple3.reactions) 
    self.assertTrue(len(mesgraph3.type_three_errors)>0) 
    #
    mesgraph4 = MESGraph(self.simple4)
    mesgraph4.analyze(self.simple4.reactions) 
    self.assertTrue(len(mesgraph4.type_four_errors)>0)  
    #
    # mesgraph5 = MESGraph(self.simple5)
    # mesgraph5.analyze(self.simple5.reactions) 
    # dih = mesgraph5.getNode(self.simple5.getMolecule(DIH))
    # din = mesgraph5.getNode(self.simple5.getMolecule(DIN))
    # self.assertTrue(mesgraph5.has_edge(dih, din))  
    # self.assertTrue(mesgraph5.has_edge(din, dih)) 
    # self.assertEqual(len(mesgraph5.type_five_errors), 1)
    #   
    mesgraph6 = MESGraph(self.simple6)
    mesgraph6.analyze(self.simple6.reactions) 
    total_errors = len(mesgraph6.type_one_errors) + \
                   len(mesgraph6.type_two_errors) + \
                   len(mesgraph6.type_three_errors) + \
                   len(mesgraph6.type_four_errors) +\
                   len(mesgraph6.type_five_errors) 
    self.assertTrue(total_errors==0)

if __name__ == '__main__':
  unittest.main()

