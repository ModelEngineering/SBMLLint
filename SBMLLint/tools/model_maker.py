"""
Creates a full Antimony model from reaction strings.
This is done by:
1. Parsing the reactions for symbols used.
2. Creating assignment statements for symbols.
3. Append these assignments to the list of reactions.
"""

import sys

EXTRANEOUS = ["+", "->", "-", "(", ")", "*", "\n", ";", ","]

class ModelMaker(object):
  """
  Creates an antimony models
  """

  def __init__(self, in_arg):
    """
    :param list-str/str reaction_strs/path:
      list-str - list of reactions
      str - file path or string of reactions
    """
    self._reaction_strs = self._getReactionstrs(in_arg)

  def _getReactionstrs(self, in_arg):
    if isinstance(in_arg, list):
      return in_arg  # list of reactions
    else:
      if "->" in in_arg:
        return in_arg.split("\n")
      else:
        with open(in_arg, "r") as fd:
          return [l for l in fd.readlines()]

  def makeModelStr(self):
    """
    Creates creates a full antimony model.
    :return str antimony_model: string of antimony model
    """
    model_strs = list(self._reaction_strs)
    [model_strs.append("%s = 0" % s) for s in self.extractSymbols()]
    return '\n'.join(model_strs)
 
  def makeModelFile(self, out_path):
    """
    Creates creates a full antimony model.
    :param str out_path: File that will be written with the model
    """
    modelstr = self.makeModelStr()
    with open(out_path, "w") as fd:
      fd.write(modelstr)

  def extractSymbols(self):
    """
    Extracts symbols from an antimony strings.
    :param list-str reaction_strs: List of reaction strings
    :return list-str: list of unique symbols found
    """
    symbols = []
    #
    for line in self._reaction_strs:
      # Remove the reaction label
      pieces = line.split(":")
      if len(pieces) == 2:
        piece = pieces[1]
      else:
        piece = pieces[0]
      for stg in EXTRANEOUS:
        piece = piece.replace(stg, " ")
      line_symbols = piece.split(" ")
      line_symbols = [s for s in line_symbols if len(s) > 0]
      for sym in line_symbols:
        try:
          float(sym)
        except:
          symbols.append(sym)
    # Get unique symbols and sort
    symbols = list(set(symbols))
    symbols.sort()
    return symbols


if __name__ == '__main__':
  in_arg = sys.stdin.read()
  maker = ModelMaker(in_arg)
  sys.stdout.write(maker.makeModelStr())
