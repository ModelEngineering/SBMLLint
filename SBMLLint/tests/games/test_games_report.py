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
G2R_CREATION = "G2R_Creation"
G2K = "G2K"
G2R = "G2R"
PG2 = "PG2"
PG2R = "PG2R"
# For Type II Error
MTR = "MTR"
MTHFR = "MTHFR"
CH3FH4 = "CH3FH4"
FH4 = "FH4"
# For Echelon, Type III Error
STATPHOSPHRYLATION = "statPhosphorylation"
PSTATDIMERISATION = "PstatDimerisation"
PSTATDIMERISATIONNUC = "PstatDimerisationNuc"


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
    self.simple3 = SimpleSBML()
    self.simple4 = SimpleSBML()
    # BIOMD0000000248 - canceling_error
    self.simple1.initialize(cn.TEST_FILE_GAMESREPORT1)
    # BIOMD0000000007 - Type I error
    self.simple2.initialize(cn.TEST_FILE_GAMESREPORT2)
    # BIOMD0000000018 - Type II error
    self.simple3.initialize(cn.TEST_FILE_GAMES_PP2)
    # BIOMD0000000167 - Echelon, Type III error
    self.simple4.initialize(cn.TEST_FILE_GAMESREPORT3)

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

  def testGetMoleculeInequalityPathReport(self):
  	if IGNORE_TEST:
  	  return
  	m = GAMES_PP(self.simple2)
  	m.analyze(error_details=False)
  	gr = GAMESReport(m)
  	count, report1 = gr.getMoleculeInequalityPathReport(G2K, PG2R, ["G2R_Creation"], 0, explain_details=False)
  	self.assertEqual(count, 1)
  	self.assertEqual(report1,
  		             "1. G2R_Creation: G2K + R -> G2R\n")
  	count, report2 = gr.getMoleculeInequalityPathReport(G2K, PG2R, ["G2R_Creation"], 0, explain_details=True)
  	self.assertEqual(count, 1)
  	self.assertEqual(report2,
  		             "G2K < PG2R by reaction(s):\n1. G2R_Creation: G2K + R -> G2R\n")

  def testReportTypeOneError(self):
  	if IGNORE_TEST:
  	  return
  	m = GAMES_PP(self.simple2)
  	m.analyze(error_details=False)
  	gr = GAMESReport(m)
  	error = [cn.PathComponents(node1=G2K, node2=G2R, reactions=[G2R_CREATION])]
  	report, error_num = gr.reportTypeOneError(error, explain_details=True)
  	self.assertEqual(error_num, [2])
  	extended_report = NULL_STR
  	extended_report = extended_report + "We detected a mass imbalance from the following reactions:\n\n"
  	extended_report = extended_report + "\nG2K = G2R by reaction(s):\n1. Rum1DegInG2R: G2R -> G2K\n\n"
  	extended_report = extended_report + "However, G2K < G2R by reaction(s):\n2. G2R_Creation: G2K + R -> G2R\n\n"
  	extended_report = extended_report + "\n----------------------------------------------------------------------\n\n"
  	extended_report = extended_report + "\n\n**********************************************************************\n\n"
  	self.assertEqual(report, extended_report)

  def testReportTypeTwoError(self):
  	if IGNORE_TEST:
  	  return
  	m = GAMES_PP(self.simple3)
  	m.analyze(error_details=False)
  	gr = GAMESReport(m)
  	som1 = m.getNode(self.simple3.getMolecule(CH3FH4))
  	som2 = m.getNode(self.simple3.getMolecule(FH4))
  	error = [[som1, som2]]
  	report, error_num = gr.reportTypeTwoError(error, explain_details=True)
  	self.assertEqual(error_num, [5])
  	LOC_START = 241
  	extended_report = NULL_STR
  	extended_report = extended_report + "{CH3FH4} < {CH2FH4=FFH2=FH2f=FH4} < {CH3FH4}\n\n"
  	extended_report = extended_report + "This indicates a mass conflict between reactions.\n"
  	extended_report = extended_report + "%s%s\n" % (PARAGRAPH_DIVIDER, PARAGRAPH_DIVIDER)
  	self.assertEqual(report[-LOC_START:], extended_report)

  def testConvertOperationSeriesToReactionOperations(self):
  	if IGNORE_TEST:
  	  return
  	m = GAMES_PP(self.simple4)
  	m.analyze(error_details=False)
  	gr = GAMESReport(m)
  	op = pd.Series([1.0, 0.5, 0.0], index = [STATPHOSPHRYLATION, PSTATDIMERISATION, PSTATDIMERISATIONNUC])
  	ro = gr.convertOperationSeriesToReactionOperations(op)
  	self.assertEqual(len(ro), 2)
  	self.assertEqual(ro[0].reaction, STATPHOSPHRYLATION)
  	self.assertEqual(ro[0].operation, 1.0)
  	self.assertEqual(ro[1].reaction, PSTATDIMERISATION)
  	self.assertEqual(ro[1].operation, 0.5)
  	# self.assertEqual(ro[2].reaction, PSTATDIMERISATIONNUC)
  	# self.assertEqual(ro[2].operation, 1.0)






