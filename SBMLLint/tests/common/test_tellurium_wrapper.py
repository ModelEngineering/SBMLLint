from SBMLLint.common import constants as cn
from SBMLLint.common.tellurium_wrapper import TelluriumWrapper


import numpy as np
import os
import unittest


IGNORE_TEST = False
INPUT = """
Hello World!
Hello World2!
"""
NUM_S1 = 2
NUM_S2 = 3
IGNORE_TEST = False
ANTIMONY_STG = '''
%dS1 -> %dS2; 1
S1 = 0
S2 = 0
''' % (NUM_S1, NUM_S2)


#############################
# Tests
#############################
class TestTelluriumWrapper(unittest.TestCase):

  def testRun(self):
    wrapper = TelluriumWrapper()
    wrapper.run("echo", INPUT)
    self.assertEqual(wrapper.return_code, 0)
    self.assertEqual(wrapper.output, INPUT)

  def testGetSBMLFromAntimony(self):
    wrapper = TelluriumWrapper()
    wrapper.run("getSBMLFromAntimony", ANTIMONY_STG)
    self.assertEqual(wrapper.return_code, 0)
    import pdb; pdb.set_trace()


if __name__ == '__main__':
  unittest.main()
