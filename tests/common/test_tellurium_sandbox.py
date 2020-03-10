from SBMLLint.common import constants as cn
from SBMLLint.common.tellurium_sandbox import TelluriumSandbox
from SBMLLint.common import tellurium_sandbox as ts
from SBMLLint.common import exceptions


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
class TestTelluriumSandbox(unittest.TestCase):

  def testRun(self):
    try:
      sandbox = TelluriumSandbox()
      sandbox.run("echo", INPUT)
      self.assertEqual(sandbox.return_code, 0)
      self.assertEqual(sandbox.output, INPUT)
    except exceptions.MissingTelluriumError:
      pass

  def testInstalledPackages(self):
    pkgs = ts.getInstalledPackages()
    self.assertTrue("scipy" in pkgs)
  
  def testGetSBMLFromAntimony(self):
    try:
      sandbox = TelluriumSandbox()
      sandbox.run("getSBMLFromAntimony", ANTIMONY_STG)
      self.assertEqual(sandbox.return_code, 0)
    except exceptions.MissingTelluriumError:
      pass
  
  def testConstructorMissingPakcage(self):
    with self.assertRaises(exceptions.MissingTelluriumError):
      sandbox = TelluriumSandbox(dependencies=["dummy"])


if __name__ == '__main__':
  unittest.main()
