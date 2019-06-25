"""
Test for Games Report
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule, MoleculeStoichiometry
from SBMLLint.common.reaction import Reaction
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.games.games_pp import SOMStoichiometry, SOMReaction, GAMES_PP
from SBMLLint.games.games_report import GAMESReport, SimplifiedReaction, NULL_STR, NUM_STAR, PARAGRAPH_DIVIDER, REPORT_DIVIDER
from SBMLLint.games.som import SOM
from SBMLLint.common import simple_sbml

import numpy as np
import pandas as pd
import os
import re
import tesbml
import unittest

IGNORE_TEST = False
# For SimplifiedReaction
CREATINEKINASE = "CreatineKinase"
PCR = "PCr"
CR = "Cr"
# For Type I Error
CDC2PHOS = "Cdc2Phos"
RUM1DEGINPG2R = "Rum1DegInPG2R"
G2K = "G2K"
PG2 = "PG2"
PG2R = "PG2R"


#############################
# Tests
#############################
class TestSimplifiedReaction(unittest.TestCase):

  def setUp(self):
  	self.simple = SimpleSBML()
  	# BIOMD0000000248 - originally for canceling_error
  	self.simple.initialize(cn.TEST_FILE_GAMESREPORT1)
  	self.mesgraph = GAMES_PP(self.simple)
  	self.mesgraph.analyze(error_details=False)
  	# Construct SimplifiedReaction
  	self.reaction = self.simple.getReaction(CREATINEKINASE)
  	self.simplified_reaction = SimplifiedReaction(self.reaction.reactants,
  		                                          self.reaction.products,
  		                                          self.reaction.label,
  		                                          self.mesgraph)

  def testConstructor(self):
  	if IGNORE_TEST:
  	  return
  	self.assertEqual(type(self.simplified_reaction.reactants[0]), MoleculeStoichiometry)
  	self.assertEqual(type(self.simplified_reaction.products[0]), MoleculeStoichiometry)
  	self.assertEqual(self.simplified_reaction.label, CREATINEKINASE)
  	self.assertEqual(type(self.simplified_reaction.mesgraph), GAMES_PP)
  	self.assertEqual(self.simplified_reaction.makeIdentifier(),
  		             self.reaction.makeIdentifier(is_include_kinetics=False))

  def testReduceBySOMs(self):
  	if IGNORE_TEST:
  	  return
  	self.simplified_reaction.reduceBySOMs()
  	self.assertEqual(len(self.simplified_reaction.reactants), 1)
  	self.assertEqual(len(self.simplified_reaction.products), 1)
  	self.assertTrue(self.simplified_reaction.reactants[0].molecule.name == PCR)
  	self.assertTrue(self.simplified_reaction.products[0].molecule.name == CR)


class TestGAMESReport(unittest.TestCase):

  def setUp(self):
    self.simple1 = SimpleSBML()
    self.simple2 = SimpleSBML()
    # BIOMD0000000248 - canceling_error
    self.simple1.initialize(cn.TEST_FILE_GAMESREPORT1)
    # BIOMD0000000007 - Type I error
    self.simple2.initialize(cn.TEST_FILE_GAMESREPORT2)

  def testReportCancelingError(self):
  	if IGNORE_TEST:
  	  return
  	m = GAMES_PP(self.simple1)
  	m.analyze(error_details=False)
  	gr = GAMESReport(m)
  	report, error_num = gr.reportCancelingError(m.canceling_errors, explain_details=True)
  	extended_report = NULL_STR
  	extended_report = extended_report + "We detected a mass imbalance from the following reactions:\n\n"
  	extended_report = extended_report + "1. OxidativePhosphorylation: 6.00 ADP + CTtis -> 6.00 ATP\n\n"
  	extended_report = extended_report + "2. ATPase: ATP -> ADP\n\n"
  	extended_report = extended_report + "*ATP and ADP have the same mass according to the above reaction\n\n\n"
  	extended_report = extended_report + "Therefore, they will result in empty product with zero mass:\n\n"
  	extended_report = extended_report + "OxidativePhosphorylation: CTtis -> \n\n"
  	extended_report = extended_report + "This indicates a mass conflict between reactions."
  	extended_report = extended_report + "\n%s%s\n" % (PARAGRAPH_DIVIDER, PARAGRAPH_DIVIDER)
  	extended_report = extended_report + "\n%s\n" % REPORT_DIVIDER
  	print("len auto error", len(report))
  	print("len manual error", len(extended_report))
  	for idx, val in enumerate(report):
  	  if val!=extended_report[idx]:
  	    print("index: ", idx)
  	    print("automatic report says::::", report[idx-20:idx+20])
  	    print("manual report says::::", extended_report[idx-20:idx+20])
  	    break
  	self.assertEqual(extended_report, report)
  	self.assertEqual(error_num, [2])

  def testGetMoleculeEqualityPath(self):
    if IGNORE_TEST:
      return
    m = GAMES_PP(self.simple2)
    m.analyze(error_details=False)
    gr = GAMESReport(m)
    som = m.getNode(self.simple2.getMolecule(G2K))
    equality_path = gr.getMoleculeEqualityPath(som, G2K, PG2R)
    self.assertTrue(len(equality_path) == 2)
    self.assertEqual(type(equality_path[0]), cn.PathComponents)
    self.assertEqual(equality_path[0].node1, G2K)
    self.assertEqual(equality_path[0].reactions, [CDC2PHOS])
    self.assertEqual(equality_path[1].node2, PG2R)
    self.assertEqual(equality_path[1].reactions, [RUM1DEGINPG2R]) 

  def testGetMoleculeEqualityPathReport(self):
    if IGNORE_TEST:
      return
    m = GAMES_PP(self.simple2)
    m.analyze(error_details=False)
    gr = GAMESReport(m)
    count, report1 = gr.getMoleculeEqualityPathReport(G2K, PG2R, 0, explain_details=False)
    self.assertEqual(count, 2)
    self.assertEqual(report1,
    	             "1. Cdc2Phos: G2K -> PG2\n2. Rum1DegInPG2R: PG2R -> PG2\n")
    count, report2 = gr.getMoleculeEqualityPathReport(G2K, PG2R, 0, explain_details=True)
    self.assertEqual(count, 2)
    self.assertEqual(report2,
    	             "\nG2K = PG2 by reaction(s):\n1. Cdc2Phos: G2K -> PG2\n\nPG2 = PG2R by reaction(s):\n2. Rum1DegInPG2R: PG2R -> PG2\n")







