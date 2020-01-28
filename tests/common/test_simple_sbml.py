"""
Tests for simple_sbml
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common import simple_sbml
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common import util

import numpy as np
import os
import libsbml
import unittest
import zipfile


IGNORE_TEST = False
NO_NAME = "dummy"
ANTIMONY_STG = '''
2S1 -> 3S2; 1
S1 = 0
S2 = 0
'''


#############################
# Tests
#############################
class TestSimpleSBML(unittest.TestCase):

  def setUp(self):
    self.simple = SimpleSBML()
    self.simple.initialize(cn.TEST_FILE)

  def testInitialize(self):
    if IGNORE_TEST:
      return
    simple = SimpleSBML()
    simple.initialize(cn.TEST_FILE)
    self.assertEqual(len(simple.reactions), cn.NUM_REACTIONS)
    self.assertEqual(len(simple.molecules), len(simple.moietys))

  def testGetMolecule(self):
    if IGNORE_TEST:
      return
    molecule1 = self.simple.molecules[0]
    name = molecule1.name
    molecule2 = self.simple.getMolecule(name)
    self.assertTrue(molecule1.isEqual(molecule2))
    self.assertIsNone(self.simple.getMolecule(NO_NAME))

  def testAdd(self):
    if IGNORE_TEST:
      return
    def test():
      self.assertEqual(len(self.simple.reactions), 2)
      for reaction in [reaction0, reaction1]:
        self.assertTrue(reaction in self.simple.reactions)
    #
    reaction0 = self.simple.reactions[0]
    reaction1 = self.simple.reactions[1]
    self.simple.reactions = [reaction0]
    self.simple.add(reaction1)
    test()
    self.simple.add(reaction1)
    test()

  def testRemove(self):
    num_reactions = len(self.simple.reactions)
    reaction0 = self.simple.reactions[0]
    reaction1 = self.simple.reactions[1]
    self.simple.remove(reaction0)
    self.assertTrue(reaction0 not in self.simple.reactions)
    self.assertTrue(reaction1 in self.simple.reactions)
    self.simple.add(reaction0)
    self.assertTrue(len(self.simple.reactions), num_reactions)

  def testGetReaction(self):
    if IGNORE_TEST:
      return
    reaction = self.simple.reactions[0]
    label = reaction.label
    reaction1 = self.simple.getReaction(label)
    self.assertTrue(reaction.isEqual(reaction1))


class TestFunctions(unittest.TestCase):

  def testReadURL(self):
    pass

  def _testIterator(self, itr):
     for item in itr:
       model = item.model
       self.assertTrue(isinstance(model.getSpecies(0),
           libsbml.Species))
     COUNT = 20
     itr = simple_sbml.modelIterator(final=COUNT)
     item_number = -1
     for item in itr:
       self.assertTrue(isinstance(item.filename, str))
       self.assertTrue(util.isSBMLModel(item.model))
       item_number = item.number
     self.assertEqual(item_number, COUNT - 1)

  def testModelIterator1(self):
    if IGNORE_TEST:
      return
    self._testIterator(simple_sbml.modelIterator(final=1))

  def testModelIterator2(self):
    if IGNORE_TEST:
      return
    self._testIterator(simple_sbml.modelIterator(
        final=1, zip_filename=None))

  def testGetZipfilePath(self):
    if IGNORE_TEST:
      return
    ffiles, zipper = simple_sbml.getZipfilePaths()
    for ffile in ffiles:
      try:
        fid = zipper.open(ffile)
        fid.close()
      except:
        assertTrue(False)
    

if __name__ == '__main__':
  unittest.main()
