from SBMLLint.common import constants as cn
from SBMLLint.tools import model_maker

import os
import unittest


IGNORE_TEST = False
TEST_IN_FILE = "test_model_maker.txt"
TEST_IN_FILE2 = "test_model_maker2.txt"
TEST_IN_PATH = os.path.join(cn.TEST_DIR, TEST_IN_FILE)
TEST_IN_PATH2 = os.path.join(cn.TEST_DIR, TEST_IN_FILE2)
TEST_OUT_FILE = "test_model_maker.ant"
TEST_OUT_PATH = os.path.join(cn.TEST_DIR, TEST_OUT_FILE)
SYM1 = "sym1"
SYM2 = "sym2"
SYM3 = "sym3"
SYM4 = "sym4"
REACTIONSTR = "%s + %s -> %s + %s" % (SYM1, SYM2, SYM3, SYM4)


#############################
# Tests
#############################
class TestModelMaker(unittest.TestCase):

  def tearDown(self):
    if os.path.isfile(TEST_OUT_PATH):
      os.remove(TEST_OUT_PATH)

  def setUp(self):
    self.maker = model_maker.ModelMaker(REACTIONSTR)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(self.maker._reaction_strs, [REACTIONSTR])
    self.maker = model_maker.ModelMaker([REACTIONSTR])
    self.assertEqual(self.maker._reaction_strs, [REACTIONSTR])
    self.maker = model_maker.ModelMaker(TEST_IN_PATH)
    self.assertGreater(len(self.maker._reaction_strs), 0)

  def testExtractSymbols(self):
    if IGNORE_TEST:
      return
    symbols = self.maker.extractSymbols()
    b = set(symbols).symmetric_difference([SYM1, SYM2, SYM3, SYM4])
    self.assertEqual(len(b), 0)
    self.maker = model_maker.ModelMaker([REACTIONSTR, REACTIONSTR])
    b = set(symbols).symmetric_difference([SYM1, SYM2, SYM3, SYM4])
    self.assertEqual(len(b), 0)

  def testMakeModelStr(self):
    if IGNORE_TEST:
      return
    model_str = self.maker.makeModelStr()
    self.assertGreater(len(model_str), len(REACTIONSTR))
   
  def testMakeModelFile(self): 
    if IGNORE_TEST:
      return
    self.assertFalse(os.path.isfile(TEST_OUT_PATH))
    model_str = self.maker.makeModelFile(TEST_OUT_PATH)
    self.assertTrue(os.path.isfile(TEST_OUT_PATH))

  def testSplitNumber(self):
    if IGNORE_TEST:
      return
    stg = "ab"
    a_str, a_int = self.maker.__class__._splitNumber(stg)
    self.assertEqual(a_str, stg)
    self.assertIsNone(a_int)
    #
    prefix = "ab"
    num = 1
    stg = "%s%d" % (prefix, num)
    a_str, a_int = self.maker.__class__._splitNumber(stg)
    self.assertEqual(a_str, prefix)
    self.assertEqual(a_int, num)
    #
    prefix = "ab"
    num = 12
    stg = "%s%d" % (prefix, num)
    a_str, a_int = self.maker.__class__._splitNumber(stg)
    self.assertEqual(a_str, prefix)
    self.assertEqual(a_int, num)
    #
    prefix = "a1b"
    num = 12
    stg = "%s%d" % (prefix, num)
    a_str, a_int = self.maker.__class__._splitNumber(stg)
    self.assertEqual(a_str, prefix)
    self.assertEqual(a_int, num)

  def testRepetitionName(self):
    if IGNORE_TEST:
      return
    def test(rname, expected, exclude_funcs=None):
      result = self.maker.__class__._makeRepetitionNames(rname,
          exclude_funcs=exclude_funcs)
      self.assertEqual(result, expected)
    #
    test("ab__bB12", "ab__bB12", exclude_funcs=[lambda n: n[0]=="b"])
    test("ab12", "ab_12")
    test("a__bb12", "a__bb_12")
    test("ab__bB12", "ab__bB_12")
    test("ab__bB12", "ab__bB_12", exclude_funcs=[lambda n: n[0]=="a"])
    test("E3SUB__SUB__misfolded__Ub2__UCHL1",
        "E3SUB__SUB__misfolded__Ub_2__UCHL1",
        exclude_funcs=[lambda n: "UCHL" in n])

  def testGetCandidateRenames(self):
    if IGNORE_TEST:
      return
    candidates = ["a_b12", "a1"]
    non_candidates = ["a1b", "a23c"]
    self.maker.symbols = list(candidates)
    self.maker.symbols.extend(non_candidates)
    result = self.maker.getCandidateRenames()
    self.assertEqual(len(candidates), len(result))

  def testGetCandidateRenames(self):
    if IGNORE_TEST:
      return
    maker = model_maker.ModelMaker(TEST_IN_PATH2)
    rename_dict1 = maker.getCandidateRenames()
    exclude_funcs = [
        lambda n: n[0]=="k",
        lambda n: n=="E1",
        lambda n: n=="E2",
        lambda n: n=="E3",
        lambda n: "UCHL1" in n,
        lambda n: "Uchl1" in n,
        lambda n: n[-2:] == "E3",
    ]
    rename_dict2 = maker.getCandidateRenames(
        exclude_funcs=exclude_funcs)
    self.assertGreater(len(rename_dict1), len(rename_dict2))

  def testReplaceSymbols(self):
    if IGNORE_TEST:
      return
    model_str = self.maker.makeModelStr()
    self.maker.replaceSymbols({SYM1: SYM2, SYM3: SYM4})
    for sym in [SYM2, SYM4]:
      self.assertEqual(self.maker.model_str.count(sym), 4)
    for sym in [SYM1, SYM3]:
      self.assertEqual(self.maker.model_str.count(sym), 0)


if __name__ == '__main__':
  unittest.main()
