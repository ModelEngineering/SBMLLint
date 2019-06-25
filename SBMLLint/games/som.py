"""Set of Molecules(SOM)."""

from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule, MoleculeStoichiometry
from SBMLLint.common.reaction import Reaction
from SBMLLint.common.simple_sbml import SimpleSBML

from collections import deque

BRACKET_OPEN = "{"
BRACKET_CLOSE = "}"

class SOM(object):
  """
  The SOM (Set Of Molecules) class represents a collection of molecules
  each of which has equal weight. Consequently, the whole molecule space 
  can be partitioned into a collection of SOMs. The uni-uni reactions
  merge multiple SOM instances to create a larger one. 
  """
  def __init__(self, molecules, reactions=set()):
    """
    :param set-Molecule molecules:
    :param set-Reaction reactions:
    """
    self.molecules = molecules
    self.reactions = reactions
    self.identifier = self.makeId()

  def __repr__(self):
    return self.identifier        
      
  def makeId(self):
    """
    Creates an identifier for the SOM to uniquely
    identifies its elements.
    :return str:
    """
    def joinMoleculeNames(molecules):
      names = [m.name for m in molecules]
      names.sort()
      return '='.join(names)
    #
    identifier = "%s%s%s" % (
        BRACKET_OPEN, 
        joinMoleculeNames(list(self.molecules)),
        BRACKET_CLOSE
        )
    return identifier

  def merge(self, som):
    """
    Creates and returns a new SOM instance
    that has the union of molecules and reactions
    :param SOM som:
    :return SOM:
    """
    new_molecules = self.molecules.union(som.molecules)
    new_reactions = self.reactions.union(som.reactions)
    new_som = SOM(molecules=new_molecules, reactions=new_reactions)
    return new_som

