"""Checks for static errors in a model."""

from SBMLLint.common import constants as cn
from SBMLLint.tools import util
from SBMLLint.structured_names.moiety_comparator import MoietyComparator


def lint(model_reference,
    mass_balance_check="structured_names"):
  """
  Reports on errors found in a model
  :param str model_reference: file, antimony, xml
  :param str mass_balance_check: how check for mass balance
  """
  util.initialize(model_reference)
  import pdb; pdb.set_trace()
  print(MoietyComparator.analyzeReactions())


if __name__ == '__main__':
  lint(cn.TEST_FILE2)
