""" Runs an python module in a separate process.
  Usage:
    runner = Runner(module_path)
    runner.execute(list-of-args, stdin-string)
    if runner.return_code == 0:
      result = runner.output  # strings written to stdout
"""


import os.path
import subprocess


class Runner(object):
  """Runs a python module as a program."""
  
  def __init__(self, python_module):
    self._python_module = python_module
    self.return_code = None  # Return code from running process
    self.output= None  # Lines written to stdout

  def execute(self, args, stdin_str):
    """
    Runs the child process from the parent process.
    Runs in the process of the caller.
    :param list-str args: list of arguments passed
       may be []
    :param str stdin_str: input to the program passed to stdin
    """
    arg_list = ['python', self._python_module]
    arg_list.extend(args)
    process = subprocess.run(
        arg_list,
        stdout=subprocess.PIPE,
        input=stdin_str, encoding='ascii')
    self.return_code = process.returncode
    self.output = process.stdout
