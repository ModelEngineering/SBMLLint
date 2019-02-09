from SBMLLint.common import constants as cn
from SBMLLint.common.reaction import Reaction
from SBMLLint.tools import util
import tesbml


import numpy as np
import os
import unittest


IGNORE_TEST = False
NUM_S1 = 2
NUM_S2 = 3
MOLECULE1 = "S1"
MOLECULE2 = "S2"
ANTIMONY_STG = '''
%d%s-> %d%s; 1
S1 = 0
S2 = 0
''' % (NUM_S1, MOLECULE1, NUM_S2, MOLECULE2)


#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

  def testInitialize(self):
    util.initialize(ANTIMONY_STG)
    reaction = Reaction.reactions[0]
    for molecule in [MOLECULE1, MOLECULE2]:
      self.assertTrue(molecule in str(reaction))


if __name__ == '__main__':
  unittest.main()
