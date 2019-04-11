"""Mass Equality Set Structure Analysis with Gaussian Elimination (MESSAGE)"""

from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import MoleculeStoichiometry
from SBMLLint.common.reaction import Reaction
from SBMLLint.common import util
from SBMLLint.games.som import SOM
from SBMLLint.common.simple_sbml import SimpleSBML

import collections
import itertools
import networkx as nx
import numpy as np
import pandas as pd
#import scipy
from scipy.linalg import lu, inv

# BIOMD 383 will test the validity of Gaussian Elimination to find errors, 
# before proceeding to building SOMGraph


class SOMStoichiometry(object):

  def __init__(self, som, stoichiometry):
    if not isinstance(som, SOM):
      raise ValueError("First argument must be a SOM.")
    if not util.isFloat(stoichiometry):
      raise ValueError("Second argument must be a float.")
    self.som = som
    self.stoichiometry = stoichiometry

  def __repr__(self):
    return "%s * % 2.2f" % (str(self.som), self.stoichiometry)


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
    self.multimulti_reactions = []
    self.permuted_matrix = None
    # L matrix from LU decomposition (not always invertible)
    self.lower_inverse  = None
    self.rref_df = None
    # Components for SOMGraph
    super(Message, self).__init__()
    self.soms = self.initializeSOMs(self.molecules)
    # networkx method
    self.add_nodes_from(self.soms)
    self.identifier = self.makeId()
    # storing errors
    self.rref_errors = []
    self.type_one_errors = []
    self.type_two_errors = []
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
    perm, lower, upper = lu(mat_t)
    # the following is applicable to transposed stoichiometry matrix
    permuted_m = inv(perm).dot(mat_t)
    pivot_index = [list(k).index(1) for k in inv(perm)]
    new_idx_mat_t = [idx_mat_t[idx] for idx in pivot_index]
    # we save as; perm_df * lower_operation * mat_t = rref_df
    perm_df = pd.DataFrame(permuted_m,
        index=new_idx_mat_t,
        columns=mat_t.columns).T
    rref_df = pd.DataFrame(upper,
        index=new_idx_mat_t,
        columns=mat_t.columns).T
    lower_operation = pd.DataFrame(lower,
        columns=new_idx_mat_t)
    self.permuted_matrix = perm_df
    self.lower_inverse = lower_operation
    self.rref_df = rref_df
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
  def getNode(self, molecule):
    """
    Find a node(SOM) containing the given molecule.
    If no such SOM exists, return False
    :param Molecule molecule:
    :return SOM/False:
    """
    for som in list(self.nodes):
      for mole in som.molecules:
        if mole.name == molecule.name:
          return som
    return False
  #
  def mergeNodes(self, som1, som2, reaction):
    """
    Merge two nodes (SOMs).
    Update arcs if applicable. 
    :param SOM som1:
    :param SOM som2:
    :param Reaction reaction:
    :return SOM new_som:
    """
    new_som = som1.merge(som2)
    new_som.reactions.add(self.simple.getReaction(reaction.label))
    for som in [som1, som2]:
      for edge in list(self.in_edges(som)):
        remaining_som = edge[0]
        reaction_label = self.get_edge_data(edge[0], edge[1])[cn.REACTION]
        self.add_edge(remaining_som, new_som, reaction=reaction_label)  
      for edge in list(self.out_edges(som)):
        remaining_som = edge[1]
        reaction_label = self.get_edge_data(edge[0], edge[1])[cn.REACTION]
        self.add_edge(new_som, remaining_som, reaction=reaction_label) 
    self.remove_nodes_from([som1, som2])
    if not self.has_node(new_som):
      self.add_node(new_som)
    return new_som
  #
  def addMultiMultiReaction(self, reaction=None):
    """
    Add a multi-multi reaction to self.multimulti_reactions
    :param reaction Reaction:
    :return bool:
    """
    if reaction not in self.multimulti_reactions:
      self.multimulti_reactions.append(reaction)  
  #
  def processUniUniReaction(self, reaction):
    """
    Process a 1-1 reaction to merge nodes.
    If no need to merge, return None.
    :param Reaction reaction:
    """
    if reaction.category != cn.REACTION_1_1:
      pass
    else:
      reactant_som = self.getNode(reaction.reactants[0].molecule)
      product_som = self.getNode(reaction.products[0].molecule)
      if reactant_som == product_som:
        return None
      else:
        new_som = self.mergeNodes(reactant_som, product_som, reaction)
        self.identifier = self.makeId()
        return new_som
  #
  def addArc(self, arc_source, arc_destination, reaction):
    """
    Add a single arc (edge) using two SOMs and reaction.
    :param SOM arc_source:
    :param SOM arc_destination:
    :param Reaction reaction:
    """
    # if there is already a preious reaction,
    if self.has_edge(arc_source, arc_destination):
      reaction_label = self.get_edge_data(arc_source, arc_destination)[cn.REACTION]
      # if reaction.label is not already included in the attribute,
      if reaction.label not in set(reaction_label):
        reaction_label = reaction_label + [reaction.label]
    else:
      reaction_label = [reaction.label]
    # overwrite the edge with new reactions set
    self.add_edge(arc_source, arc_destination, reaction=reaction_label)
  # 
  def processUniMultiReaction(self, reaction):
    """
    Process a 1-n reaction to add arcs.
    Since the mass of reactant is greater than
    that of each product, it adds arcs by
    addArc(source=products, destination=reactant).
    :param Reaction reaction:
    """
    if reaction.category != cn.REACTION_1_n:
      pass
    else:
      destination = [reaction.reactants[0].molecule]
      source = [product.molecule for product in reaction.products]
      arcs = itertools.product(source, destination)
      for arc in arcs:
        if not self.checkTypeOneError(arc, reaction):
          som_source = self.getNode(arc[0])
          som_destination = self.getNode(arc[1])
          self.addArc(som_source, som_destination, reaction)
      self.identifier = self.makeId()
  #
  def processMultiUniReaction(self, reaction):
    """
    Process a n-1 reaction to add arcs.
    Since the mass of product is greater than
    that of each reactant, it adds arcs by
    addArc(source=reactants, destination=product).
    :param Reaction reaction:
    """
    if reaction.category != cn.REACTION_n_1:
      pass
    else:
      destination = [reaction.products[0].molecule]
      source = [reactant.molecule for reactant in reaction.reactants]
      arcs = itertools.product(source, destination)
      for arc in arcs:
        if not self.checkTypeOneError(arc, reaction):
          som_source = self.getNode(arc[0])
          som_destination = self.getNode(arc[1])
          self.addArc(som_source, som_destination, reaction)
      self.identifier = self.makeId()
  #
  def addTypeOneError(self, mole1, mole2, reaction):
    """
    Add Type I Error components to self.type_one_errors
    All components of resulting PathComponents are str
    :param Molecule mole1:
    :param Molecule mole2:
    :param Reaction reaction:
    :return bool flag:
    """
    flag = False
    for component in self.type_one_errors:
      if (component.node1==mole1.name) and (component.node2==mole2.name):
        new_component = cn.PathComponents(node1=mole1.name, 
                                          node2=mole2.name,
                                          reactions=component.reactions+[reaction.label])
        self.type_one_errors.remove(component)
        self.type_one_errors.append(new_component)
        flag = True
        break
    if not flag:
      self.type_one_errors.append(cn.PathComponents(node1=mole1.name, 
                                                    node2=mole2.name,
                                                    reactions=[reaction.label]))
      flag = True
    return flag
  #
  def checkTypeOneError(self, arc, inequality_reaction=None):
    """
    Check Type I Error of an arc.
    If both source and destination are found
    in the same SOM, send error message and return True.
    If not, return False.
    :param tuple-Molecule arc:
    :param Reaction inequality_reaction:
    :return bool:
    """
    som1 = self.getNode(arc[0])
    som2 = self.getNode(arc[1])
    if som1 == som2:
      self.addTypeOneError(arc[0], arc[1], inequality_reaction)
      return True
    else:
      return False
  #
  def addTypeTwoError(self, cycle):
    """
    Add Type II Error components to self.type_two_errors
    which is a list of lists
    All components of resulting PathComponents are str
    :param list-SOM cycle:
    """
    # exceptionally, here PathComponents are
    # node1=[], node2=[], reactions=[] and their index
    # of each component will match. All elements within nodes
    # are in the same SOM
    error_cycle = []
    for node_idx in range(len(cycle)-1):
      som1 = cycle[node_idx]
      som2 = cycle[node_idx+1]
      som1_moles = {mole.name for mole in list(som1.molecules)}
      som2_moles = {mole.name for mole in list(som2.molecules)}
      reactions = self.get_edge_data(som1, som2)[cn.REACTION]
      # all reactions (in an edge), should create a single PathComponent
      nodes1 = []
      nodes2 = []
      reaction_labels = []
      for r in reactions:
        reaction = self.simple.getReaction(r)
        if reaction.category == cn.REACTION_n_1:
          sources = {r.molecule.name for r in reaction.reactants}
          destinations = {p.molecule.name for p in reaction.products}
        elif reaction.category == cn.REACTION_1_n:
          sources = {p.molecule.name for p in reaction.products}
          destinations = {r.molecule.name for r in reaction.reactants}
        # for any reaction that addes arcs, len(nodes2)==1
        node2 = list(destinations.intersection(som2_moles))[0]
        for node1 in list(sources.intersection(som1_moles)):
          nodes1.append(node1)
          nodes2.append(node2)
          reaction_labels.append(reaction.label)
      error_cycle.append(cn.PathComponents(node1=nodes1, 
                                           node2=nodes2,
                                           reactions=reaction_labels))
    som1 = cycle[-1]
    som2 = cycle[0]
    som1_moles = {mole.name for mole in list(som1.molecules)}
    som2_moles = {mole.name for mole in list(som2.molecules)}
    reactions = self.get_edge_data(som1, som2)[cn.REACTION]
    # all reactions (in an edge), should create a single PathComponent
    nodes1 = []
    nodes2 = []
    reaction_labels = []
    for r in reactions:
      reaction = self.simple.getReaction(r)
      if reaction.category == cn.REACTION_n_1:
        sources = {r.molecule.name for r in reaction.reactants}
        destinations = {p.molecule.name for p in reaction.products}
      elif reaction.category == cn.REACTION_1_n:
        sources = {p.molecule.name for p in reaction.products}
        destinations = {r.molecule.name for r in reaction.reactants}
        # for any reaction that addes arcs, len(nodes2)==1
      node2 = list(destinations.intersection(som2_moles))[0]
      for node1 in list(sources.intersection(som1_moles)):
        nodes1.append(node1)
        nodes2.append(node2)
        reaction_labels.append(reaction.label)
    error_cycle.append(cn.PathComponents(node1=nodes1, 
                                        node2=nodes2,
                                        reactions=reaction_labels))
    self.type_two_errors.append(error_cycle)  
  #
  def checkTypeTwoError(self):
    """
    Check Type II Error (cycles) of a MESGraph.
    If there is at least one cycle, 
    report an error message, related reactions
    and return True.
    If there is no cycle, return False. 
    :return bool:
    """
    graph = nx.DiGraph()
    graph.add_edges_from(self.edges)
    cycles = list(nx.simple_cycles(graph))
    if len(cycles) == 0:
      return False
    else:
      for cycle in cycles:
        ######self.addTypeTwoError(cycle)######
        self.type_two_errors.append(cycle)
        # if not self.type_two_error:
        #   self.type_two_error = True
      return True
  #
  def processErrorReactions(self, reaction):
    """
    Simply add error reactions to error
    :param Reaction reaction:
    :return bool:
    """
    error_reactions = [r.label for r in self.rref_errors]
    if reaction.label not in error_reactions:
      self.rref_errors.append(reaction)
    ## Here, we may need to add more info that will 
    ## help us track the operations that lead to this error 
    return True
  #
  def analyze(self, reactions=None, error_details=True):
    """
    Using the stoichiometry matrix, compute
    row reduced echelon form and create SOMGraph
    Add arcs or sending error messages using
    checkTypeOneError or checkTypeTwoError.
    :param bool error_details:
    :return bool:
    """
    if reactions is None:
      reactions = self.simple.reactions
    # Associate the reaction category with the function
    # that processes that category
    report = cn.NULL_STR
    reaction_dic = {
        cn.REACTION_1_1: self.processUniUniReaction,
        cn.REACTION_1_n: self.processUniMultiReaction,
        cn.REACTION_n_1: self.processMultiUniReaction,
        cn.REACTION_n_n: self.addMultiMultiReaction,
        }
    # Process each type of reaction
    for category in reaction_dic.keys():
      for reaction in [r for r in reactions if r.category == category]:
        func = reaction_dic[category]
        func(reaction)
    #

    # Decompose matrix and prepare reduced reactions
    rref_df = self.decomposeMatrix(self.stoichiometry_matrix)
    reduced_reactions = self.convertMatrixToReactions(self.simple, rref_df)
    self.reduced_reactions = reduced_reactions
    report = cn.NULL_STR
    reaction_dic = {
        cn.REACTION_ERROR: self.processErrorReactions,
        cn.REACTION_1_1: self.processUniUniReaction,
        cn.REACTION_1_n: self.processUniMultiReaction,
        cn.REACTION_n_1: self.processMultiUniReaction
        }
    # Process each type of reaction
    for category in reaction_dic.keys():
      for reaction in [r for r in reduced_reactions if r.category == category]:
        func = reaction_dic[category]
        func(reaction)
    #
    self.checkTypeTwoError()
    # print("We just analyzed the data...")
    # print("RREF error: ", self.rref_errors)
    # print("Type I error: ", self.type_one_errors)
    # print("Type II error: " , self.type_two_errors)
    if self.rref_errors or self.type_one_errors or self.type_two_errors:
      return True
    else:
      return False








