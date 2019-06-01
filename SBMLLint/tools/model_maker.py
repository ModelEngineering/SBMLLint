"""
Creates a full Antimony model from reaction strings.
This is done by:
1. Parsing the reactions for symbols used.
2. Creating assignment statements for symbols.
3. Append these assignments to the list of reactions.
Assumes that name segments are separated by a DOUBLEUNDERSCORE.

TODO:
1. Keep state of updates to model
2. findCandidateRenames - finds names that should be renamed
   to moiety stoichiometries (e.g. *[0-9]_. Returns a dict
   of current name (key) and proposed name change (value)
3. changeName - inputs a dict of old and new name. Changes names
   from longest to shortest.
"""

import pandas as pd
import sys

EXTRANEOUS = ["+", "->", "-", "(", ")", "*", "\n", ";", ","]
UNDERSCORE = "_"
DOUBLEUNDERSCORE = "__"
IGNORE_SYMBOLS = ["pow", "/"]

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
    self.symbols = None  # symbols found
    self.model_str = None  # full model

  def _getReactionstrs(self, in_arg):
    if isinstance(in_arg, list):
      return in_arg  # list of reactions
    else:
      if "->" in in_arg:
        return in_arg.split("\n")
      else:
        with open(in_arg, "r") as fd:
          return [l for l in fd.readlines()]

  @staticmethod
  def _splitNumber(stg):
    """
    Extracts an ending integer from a string if there is one.
    :return str, int: prefix, number
    """
    prefix = stg
    num = None
    for idx in range(len(stg)):
      try:
        _ = int(prefix[-1])
        if num is None:
          num = prefix[-1]
        else:
          num = prefix[-1] + num
        prefix = prefix[:-1]
      except:
        break
    if num is None:
      return prefix, None
    else:
      return prefix, int(num)

  @classmethod
  def _makeRepetitionNames(cls, symbol, exclude_funcs=None):
    """
    Changes a symbol name to use moiety-stoichiometry repetitions 
    if it has repetition numbers. Has arguments to selectively
    exclude symbols
    :param str symbol: candidate repetition name
    :param list-BooleanFunction exclude_funcs: returns True if name
        should be excluded for consideration as a repetition count.
    :return str:
    """
    if exclude_funcs is None:
      exclude_funcs = []
    splits = symbol.split(DOUBLEUNDERSCORE)
    is_changed = False
    for idx, stg in enumerate(splits):
      is_ignore = False
      for func in exclude_funcs:
        if func(stg):
          is_ignore = True
      if not is_ignore:
        moiety, stoich = cls._splitNumber(stg)
        if stoich is not None:
          is_changed = True
          splits[idx] = "%s%s%s" % (moiety, UNDERSCORE, stoich)
    if is_changed:
      return DOUBLEUNDERSCORE.join(splits)
    else:
      return symbol

  def getCandidateRenames(self, exclude_funcs=None):
    """
    Finds symbols that should be renamed to use repetition counts.
    :param list-BooleanFunction exclude_funcs: applied to name parts
    :return dict: key is current name; value is proposed name
    """
    if self.symbols is None:
      self.extractSymbols()
    rename_dict = {s: self.__class__._makeRepetitionNames(s, exclude_funcs=exclude_funcs)
        for s in self.symbols}
    rename_dict = {k: v for k, v in rename_dict.items() if k != v}
    return rename_dict

  def replaceSymbols(self, symbol_dict, is_sort=True):
    """
    Replaces the symbol (key) with its value.
    :param dict symbol_dict:
    :param bool is_sort: sort replacement strings by size to
        avoid substrint replacements
    Updates self.model_str
    :return str:
    """
    ser = pd.Series([len(k) for k in symbol_dict.keys()])
    ser.index = [s for s in symbol_dict.keys()]
    if is_sort:
      ser = ser.sort_values(ascending=False)
    for k in ser.index:
      self.model_str = self.model_str.replace(k, symbol_dict[k])
    return self.model_str

  def makeModelStr(self):
    """
    Creates creates a full antimony model.
    :return str antimony_model: string of antimony model
    """
    model_strs = list(self._reaction_strs)
    [model_strs.append("%s = 0" % s) for s in self.extractSymbols()]
    self.model_str = '\n'.join(model_strs)
    return self.model_str
 
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
    self.symbols = []
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
          self.symbols.append(sym)
    # Get unique symbols and sort
    self.symbols = [s for s in set(self.symbols)
        if not s in IGNORE_SYMBOLS]
    self.symbols.sort()
    return self.symbols


if __name__ == '__main__':
  in_arg = sys.stdin.read()
  maker = ModelMaker(in_arg)
  sys.stdout.write(maker.makeModelStr())
