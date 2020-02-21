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
STATPHOSPHORYLATION = "statPhosphorylation"
PSTATDIMERISATION = "PstatDimerisation"
PSTATDIMERISATIONNUC = "PstatDimerisationNuc"
SPECIES_TEST = "species_test"
PSTAT_SOL = "Pstat_sol"
PSTATDIMER_NUC = "PstatDimer_nuc"
PSTAT_NUC = "Pstat_nuc"
STAT_SOL = "stat_sol"
PSTATDIMER_SOL = "PstatDimer_sol"


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
    extended_report = extended_report + "We detected a mass imbalance\n"
    extended_report = extended_report + ": OxidativePhosphorylation: CTtis -> \n\n"
    extended_report = extended_report + "from the following isolation set:\n\n"
    extended_report = extended_report + "1. OxidativePhosphorylation: 6.00 ADP + CTtis -> 6.00 ATP\n"
    extended_report = extended_report + "2. ATPase: ATP -> ADP\n"
    extended_report = extended_report + "*ATP and ADP have the same mass according to the above reaction\n"
    extended_report = extended_report + "\n%s%s\n" % (PARAGRAPH_DIVIDER, PARAGRAPH_DIVIDER)
    extended_report = extended_report + "\n%s\n" % REPORT_DIVIDER
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
  	op = pd.Series([1.0, 0.5, 0.0], index = [STATPHOSPHORYLATION, PSTATDIMERISATION, PSTATDIMERISATIONNUC])
  	ro = gr.convertOperationSeriesToReactionOperations(op)
  	self.assertEqual(len(ro), 2)
  	self.assertEqual(ro[0].reaction, STATPHOSPHORYLATION)
  	self.assertEqual(ro[0].operation, 1.0)
  	self.assertEqual(ro[1].reaction, PSTATDIMERISATION)
  	self.assertEqual(ro[1].operation, 0.5)

  def testGetOperationMatrix(self):
    if IGNORE_TEST:
      return
    m1 = GAMES_PP(self.simple1)
    m1.analyze(error_details=False)
    gr1 = GAMESReport(m1)
    self.assertTrue(gr1.getOperationMatrix() is None)
    m4 = GAMES_PP(self.simple4)
    m4.analyze(error_details=False)
    gr4 = GAMESReport(m4)
    op_mat = gr4.getOperationMatrix()
    self.assertEqual(op_mat.loc[STATPHOSPHORYLATION, STATPHOSPHORYLATION], 1.0)
    self.assertEqual(op_mat.loc[PSTATDIMERISATION, PSTATDIMERISATION], 1.0)
    self.assertEqual(op_mat.loc[PSTATDIMERISATIONNUC, PSTATDIMERISATIONNUC], 1.0)
    self.assertEqual(op_mat.loc[PSTATDIMERISATIONNUC, STATPHOSPHORYLATION], 0.0)
    self.assertEqual(op_mat.loc[STATPHOSPHORYLATION, PSTATDIMERISATIONNUC], -0.5)

  def testGetResultingSeries(self):
    if IGNORE_TEST:
      return
    m = GAMES_PP(self.simple4)
    m.analyze(error_details=False)
    gr = GAMESReport(m)
    resulting_series = gr.getResultingSeries(STATPHOSPHORYLATION)
    print(resulting_series)
    self.assertEqual(resulting_series["{" + SPECIES_TEST + "}"], 1.0)
    self.assertEqual(resulting_series["{" + PSTAT_SOL + "}"], 0.0)
    self.assertEqual(resulting_series[m.getNode(m.simple.getMolecule(PSTATDIMER_NUC)).identifier], 0.0)
    self.assertEqual(resulting_series[m.getNode(m.simple.getMolecule(PSTAT_NUC)).identifier], 0.0)

  def testGetOperationStoichiometryMatrix(self):
    if IGNORE_TEST:
      return
    m = GAMES_PP(self.simple4)
    m.analyze(error_details=False)
    gr = GAMESReport(m)
    op = pd.Series([1.0, 0.5, 0.0], index = [STATPHOSPHORYLATION, PSTATDIMERISATION, PSTATDIMERISATIONNUC])
    ro = gr.convertOperationSeriesToReactionOperations(op)
    osm = gr.getOperationStoichiometryMatrix(ro)
    self.assertEqual(osm.loc[SPECIES_TEST, STATPHOSPHORYLATION], 1.0)
    self.assertEqual(osm.loc[SPECIES_TEST, PSTATDIMERISATION], 0.0)
    self.assertEqual(osm.loc[PSTAT_SOL, STATPHOSPHORYLATION], 1.0)
    self.assertEqual(osm.loc[PSTAT_SOL, PSTATDIMERISATION], -2.0)

  def testGeInferredReaction(self):
    if IGNORE_TEST:
      return
    m = GAMES_PP(self.simple4)
    m.analyze(error_details=False)
    gr = GAMESReport(m)
    op = pd.Series([1.0, 0.5, 0.0], index = [STATPHOSPHORYLATION, PSTATDIMERISATION, PSTATDIMERISATIONNUC])
    ro = gr.convertOperationSeriesToReactionOperations(op)
    inferred_reaction = gr.getInferredReaction(ro)
    self.assertEqual(len(inferred_reaction.reactants), 1)
    self.assertEqual(len(inferred_reaction.products), 2) 
    self.assertEqual(inferred_reaction.reactants[0].molecule.name, STAT_SOL)
    self.assertTrue(inferred_reaction.products[0].molecule.name in {PSTATDIMER_SOL, SPECIES_TEST})
    self.assertTrue(inferred_reaction.products[1].molecule.name in {PSTATDIMER_SOL, SPECIES_TEST})
    self.assertEqual(inferred_reaction.reactants[0].stoichiometry, 1.0)
    self.assertEqual({p.stoichiometry for p in inferred_reaction.products}, {0.5, 1.0})

  def testReportReactionsInSOM(self):
    if IGNORE_TEST:
      return
    m = GAMES_PP(self.simple4)
    m.analyze(error_details=False)
    gr = GAMESReport(m)
    som = m.getNode(m.simple.getMolecule(PSTATDIMER_NUC))
    report, error_num = gr.reportReactionsInSOM(som, 0)
    common_part = "1. PstatDimer__import: PstatDimer_sol -> PstatDimer_nuc\n"
    self.assertEqual(error_num, 1)
    self.assertTrue(report == common_part)

  def testReportEchelonError(self):
  	if IGNORE_TEST:
  	  return
  	m = GAMES_PP(self.simple4)
  	m.analyze(error_details=False)
  	gr = GAMESReport(m)
  	report, error_num = gr.reportEchelonError(m.echelon_errors, explain_details=True)
  	self.assertEqual(error_num, [3])
  	extended_report = NULL_STR
  	extended_report = extended_report + "will result in empty reactant with zero mass:\n\n:  -> {species_test}\n\n"
  	extended_report = extended_report + "\n----------------------------------------------------------------------\n"
  	extended_report = extended_report + "\n----------------------------------------------------------------------\n\n"
  	extended_report = extended_report + "\n\n**********************************************************************\n\n"
  	self.assertEqual(report[-288:], extended_report)

  def testReportTypeThreeError(self):
    if IGNORE_TEST:
      return
    m = GAMES_PP(self.simple4)
    m.analyze(error_details=False)
    gr = GAMESReport(m)
    report, error_num = gr.reportTypeThreeError(m.type_three_errors, explain_details=True)
    self.assertEqual(error_num, [3])
    pseudo_inequality_report = NULL_STR
    pseudo_inequality_report = pseudo_inequality_report + "6. statPhosphorylation: stat_sol -> Pstat_sol + species_test\n"
    pseudo_inequality_report = pseudo_inequality_report + "(pseudo 6.) statPhosphorylation: {Pstat_nuc=stat_nuc=stat_sol} -> "
    pseudo_inequality_report1 = pseudo_inequality_report + "{species_test} + {Pstat_sol}"
    pseudo_inequality_report2 = pseudo_inequality_report + "{Pstat_sol} + {species_test}"
    inference_report1 = "the masses of {Pstat_sol} and {Pstat_nuc=stat_nuc=stat_sol} are unequal."
    inference_report2 = "the masses of {Pstat_nuc=stat_nuc=stat_sol} and {Pstat_sol} are unequal."
    self.assertTrue(report[-460:-305]==pseudo_inequality_report1 or report[-460:-305]==pseudo_inequality_report2)
    self.assertTrue(report[-293:-221]==inference_report1 or report[-293:-221]==inference_report2)



