""" Wraps tellurium methods to run in anothe process.  """

import argparse
import subprocess
import sys


class TelluriumWrapper(object):
  """
  Runs python code in a separate process (sandbox).
  Note that there are two instances of this class that are
  created; one in the parent process (which uses the main and run
  methods) and one in the child (which uses the execute method).
  Usage:
    wrapper.run(method_name, string_argument)
  Outputs
    wrapper.return_code - should be 0
    wrapper.output - string result
  """

  def __init__(self):
    self.return_code = None
    self.output = None

  def echo(self, input_string):
    """
    Executes the sandboxed python codes in the child process.
    :param str input_string:
    :return str:
    Notes:
      1. Writes double newlines ("\n")
    """
    sys.stdout.writelines(input_string)

  def getSBMLFromAntimony(self, input_string):
    import tellurium as te
    #
    rr = te.loada(input_string)
    sbml = rr.getSBML()
    sys.stdout.writelines(sbml)

  def run(self, method, input_string):
    """
    Runs the child process from the parent process.
    Strips double newlines.
    """
    process = subprocess.run(['python', __file__, method],
        stdout=subprocess.PIPE,
        input=input_string, encoding='ascii')
    self.return_code = process.returncode
    self.output = process.stdout

  def main(self):
    """
    The first argument is the method to call.
    """
    cmd = "self.%s(sys.stdin)" % sys.argv[1]
    exec(cmd)


if __name__ == '__main__':
  wrapper = TelluriumWrapper()
  wrapper.main()
