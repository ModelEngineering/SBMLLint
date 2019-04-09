"""Mass Equality Set Structure Analysis with Gaussian Elimination (MESSAGE)"""

from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import MoleculeStoichiometry
from SBMLLint.common.reaction import Reaction
from SBMLLint.games.som import SOM
from SBMLLint.common.simple_sbml import SimpleSBML

import collections
import itertools
import networkx as nx
import numpy as np
import pandas as pd
#import scipy
from sympy.matrices import Matrix, eye

# BIOMD 383 will test the validity of Gaussian Elimination to find errors, 
# before proceeding to building MESGraph



class Message(nx.DiGraph):
  """
  Similar to MESGraph, The Message algorithm creates
  a directed graph of SOMs, but before creating the graph
  it creates a row-reduced echolon metrix using the 
  stoichiometry matrix. 
  Type I Error occurs when we find inequality between two molecules
  in the same SOM, because each element in a SOM has the same weight.
  Type II Error implies there is cyclism between molecules, such as
  A < B < C < ... < A, which is physically impossible.
  """
  def __init__(self, simple=None):
    """
    :param SimpleSBML simple:
    """
    self.simple = simple
    self.reactions = self._getNonBoundaryReactions(simple)
    self.molecules = self._getNonBoundaryMolecules(simple, self.reactions)
    self.stoichiometry_matrix = self.getStoichiometryMatrix(self.reactions, self.molecules)
    self.reduced_reactions = []
    # INVERSE of L from LU decomposition
    self.lower_inverse = None
    self.permuted_matrix = None
    # Components for SOMGraph
    super(Message, self).__init__()
    self.soms = self.initializeSOMs(self.molecules)
    # networkx method
    self.add_nodes_from(self.soms)
    self.identifier = self.makeId()

  #
  def __repr__(self):
    return self.identifier
  #
  def makeId(self):
    """
    Construct an identifier for the graph.
    :return str:
    """
    identifier = ""
    if self.edges:
      for edge in self.edges:
        identifier = identifier + str(edge[0]) + cn.ARC_ARROW + str(edge[1]) + "\n"
    for key, node in enumerate(nx.isolates(self)):
      identifier = identifier + str(node)
      if key < (len(list(nx.isolates(self)))-1):
          identifier = identifier + cn.KINETICS_SEPARATOR
    # Return the identifier
    return identifier
  #
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
  #
  def _getNonBoundaryMolecules(self, simple, reactions):
    """
    Get list of non-boundary molecules
    :param SimpleSBML simple:
    :return list-Molecule:
    """
    molecules = set()
    for reaction in reactions:
      reactants = {simple.getMolecule(r.molecule.name) for r in reaction.reactants}
      products = {simple.getMolecule(r.molecule.name) for r in reaction.products}
      molecules = molecules.union(reactants)
      molecules = molecules.union(products)
    return list(molecules)
  #
  def initializeSOMs(self, molecules):
    """
    Create a list of one-molecule SOMs
    :param list-Molecule molecules:
    :return list-SOM soms:
    """
    soms = []
    for molecule in molecules:
      if molecule.name == cn.EMPTYSET:
        continue
      else:
        soms.append(SOM({molecule}))
    return soms
  #
  def getStoichiometryMatrix(self, reactions, molecules):
    """
    Creates a full stoichiometry matrix
    using non-boundary reactions.
    Helped by https://gist.github.com/lukauskas/d1e30bdccc5b801d341d
    :return pd.DataFrame:
    """
    reaction_labels = [r.label for r in reactions]
    molecule_names = [m.name for m in molecules]
    stoichiometry_matrix = pd.DataFrame(0.0, index=molecule_names, columns=reaction_labels)
    for reaction in reactions:
      reactants = {r.molecule.name:r.stoichiometry for r in reaction.reactants}
      products = {p.molecule.name:p.stoichiometry for p in reaction.products}
      reaction_molecules = list(set(reactants.keys()).union(products.keys()))
      for molecule_name in reaction_molecules:
        net_stoichiometry = products.get(molecule_name, 0.0) - reactants.get(molecule_name, 0.0)
        stoichiometry_matrix[reaction.label][molecule_name] = net_stoichiometry
    return stoichiometry_matrix
  #
  def decomposeMatrix(self, mat_df):
    """
    LU decomposition of the stoichiometry matrix.
    First it transposes the input matrix
    and find L, U matrices and permutation list. 
    :param pandas.DataFrame mat_df:
    :yield pandas.DataFrame new_df:
    """
    mat_t = mat_df.T
    idx_mat_t = mat_t.index
    # LU decomposition
    m = Matrix(mat_df.T)
    lower, upper, perm = m.LUdecomposition()
    permuted_m = m.permuteFwd(perm)
    new_idx_mat_t = list(Matrix(idx_mat_t).permuteFwd(perm))
    perm_df = pd.DataFrame(np.array(permuted_m).astype(np.float64),
        index=new_idx_mat_t,
        columns=mat_t.columns).T
    rref_df = pd.DataFrame(np.array(upper).astype(np.float64),
        index=new_idx_mat_t,
        columns=mat_t.columns).T
    lower_operation = pd.DataFrame(np.array(lower.inv()).astype(np.float64),
        columns=new_idx_mat_t)
    self.permuted_matrix = perm_df
    self.lower_inverse = lower_operation
    return rref_df
  #
  def getReactionSummaryCategory(self, reactants, products):
    """
    Return category of reaction. Return reaction_n_n
    if none of the above applies
    :param list-MoleculeStoichiometry reactants:
    :param list-Moleculestoichiometry products:
    :return str reaction_category:
    """
    num_reactants = len([r.molecule for r in reactants \
                         if r.molecule.name!=cn.EMPTYSET])
    num_products = len([p.molecule for p in products \
                        if p.molecule.name!=cn.EMPTYSET])
    stoichiometry_reactants = [r.stoichiometry for r \
                                  in reactants \
                                  if r.molecule.name!=cn.EMPTYSET]
    stoichiometry_products = [p.stoichiometry for p \
                             in products \
                             if p.molecule.name!=cn.EMPTYSET]
    for reaction_category in cn.REACTION_SUMMARY_CATEGORIES:
      if reaction_category.predicate(num_reactants, num_products, 
                                     stoichiometry_reactants, 
                                     stoichiometry_products):
        return reaction_category.category
    # if none of the above, return reaction_n_n
    return cn.REACTION_n_n
  #
  def convertMatrixToReactions(self, simple, mat_df):
    """
    Convert a stoichiometry matrix, 
    where columns are reactions and 
    rows are molecules(species),
    to simpleSBML reactions. 
    :param simpleSBML simple:
    :param pandas.DataFrame mat_df:
    :return list-ReactionSummary reactions:
    """
    reactions = []
    for reaction_name in mat_df.columns:
      reaction = simple.getReaction(reaction_name)
      reduced_reaction_series = mat_df[reaction_name]
      reactants = [MoleculeStoichiometry(simple.getMolecule(molecule), 
                                     abs(reduced_reaction_series[molecule])) \
              for molecule in reduced_reaction_series.index if reduced_reaction_series[molecule]<0]
      products = [MoleculeStoichiometry(simple.getMolecule(molecule), 
                                     reduced_reaction_series[molecule]) \
              for molecule in reduced_reaction_series.index if reduced_reaction_series[molecule]>0]
      reactions.append(cn.ReactionSummary(label=reaction_name, 
                                      reactants=reactants,
                                      products=products,
                                      category=self.getReactionSummaryCategory(reactants, products)))
    return reactions
  #
  def analyze(self, error_details=True):
    """
    Using the stoichiometry matrix, compute
    row reduced echelon form and create SOMGraph
    Add arcs or sending error messages using
    checkTypeOneError or checkTypeTwoError.
    :param bool error_details:
    :return str:
    """
    report = NULL_STR
    pass








