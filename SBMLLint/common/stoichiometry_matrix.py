"""Stoichiometry Matrix."""

from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule, MoleculeStoichiometry
from SBMLLint.common.reaction import Reaction
from SBMLLint.games.som import SOM
from SBMLLint.common.simple_sbml import SimpleSBML

import numpy as np
import pandas as pd
from scipy.optimize import linprog
import warnings


class StoichiometryMatrix(object):
  """
  Creates a full stoichiometry matrix from simpleSBML 
  that correctly incorporates boundary species 
  and use linear programming to determine stoichiometric consistency.
  """    
  def __init__(self, simple=None):
    self.reactions = self._getNonBoundaryReactions(simple)
    self.molecules = self._getNonBoundaryMolecules(simple)
    self.stoichiometry_matrix = self.makeStoichiometryMatrix()
    self.consistent = None
    self.result = None

  def _getNonBoundaryReactions(self, simple):
    """
    Get list of non-boundary reacetions
    :param SimpleSBML simple:
    :return list-Reaction:
    """
    reactions = []
    for reaction in simple.reactions:
      if reaction.category != cn.REACTION_BOUNDARY:
        reactions.append(reaction)
    return reactions

  def _getNonBoundaryMolecules(self, simple):
    """
    Get list of non-boundary molecules
    :param SimpleSBML simple:
    :return list-Molecule.name:
    """
    molecules = set()
    for reaction in self.reactions:
      reactants = {r.molecule.name for r in reaction.reactants}
      products = {r.molecule.name for r in reaction.products}
      molecules = molecules.union(reactants)
      molecules = molecules.union(products)
    return list(molecules)

  def makeStoichiometryMatrix(self):
    """
    Creates a full stoichiometry matrix
    using non-boundary reactions.
    Helped by https://gist.github.com/lukauskas/d1e30bdccc5b801d341d
    :return pd.DataFrame:
    """
    reaction_labels = [r.label for r in self.reactions]
    stoichiometry_matrix = pd.DataFrame(0.0, index=self.molecules, columns=reaction_labels)
    for reaction in self.reactions:
      reactants_raw = [(r.molecule.name, r.stoichiometry) for r in reaction.reactants]
      products_raw = [(p.molecule.name, p.stoichiometry) for p in reaction.products]
      reactants_key = list(set([r[0] for r in reactants_raw]))
      products_key = list(set([p[0] for p in products_raw]))
      # Below is to deal with reactions with repeated species; e.g., S0 -> S1 + S1
      reactants = {r_k:sum([r[1] for r in reactants_raw if r[0]==r_k]) for r_k in reactants_key}
      products = {p_k:sum([p[1] for p in products_raw if p[0]==p_k]) for p_k in products_key}
      reaction_molecules = list(set(reactants.keys()).union(products.keys()))
      for molecule_name in reaction_molecules:
        net_stoichiometry = products.get(molecule_name, 0.0) - reactants.get(molecule_name, 0.0)
        stoichiometry_matrix[reaction.label][molecule_name] = net_stoichiometry
    return stoichiometry_matrix

  def isConsistent(self, is_report_warning=True):
    """
    Runs linear programmming (LP) to determine 
    stoichiometric inconsistency. 
    If consistent return True,
    else return False. 
    :param bool is_report_warning: report optimization warnings
    :return bool:
    """
    s_matrix_t = self.stoichiometry_matrix.T
    # number of reactions
    nreac = s_matrix_t.shape[0]
    # number of chemical species
    nmet = s_matrix_t.shape[1]
    #
    b = np.zeros(nreac)
    c = np.ones(nmet)
    # Linear programming. c is constraint (here, zero), 
    # b is vector of possible values for molecule vector. 
    if not is_report_warning:
      warnings.simplefilter("ignore")
    try:
      res = linprog(c, A_eq=s_matrix_t, b_eq=b, bounds=(1, None))
      self.result = res
      is_success = True
    except:
      is_success = False
    if not is_success:
      msg = "*** Failed to solve the stoichiometry matrix."
      raise RuntimeError(msg)
    if res.status == 0:
      self.consistent = True
    else:
      self.consistent = False
    #
    return self.consistent
















