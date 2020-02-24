"""GAMES Plus (Reduced) Row Echelon Form (GAMES_PP)"""

from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule, MoleculeStoichiometry
from SBMLLint.common.reaction import Reaction
from SBMLLint.common import util
from SBMLLint.games.som import SOM
from SBMLLint.common.simple_sbml import SimpleSBML

import collections
import itertools
import networkx as nx
import numpy as np
import pandas as pd
from scipy.linalg import lu, inv

GAMESErrors = collections.namedtuple("GAMESErrors", 
    "type_one")
ErrorSummary = collections.namedtuple("ErrorSummary",
    "type errors")
TOLERANCE = 0.0001
TYPE_I = "type1"
TYPE_II = "type2"
TYPE_III = "type3"
CANCELING = "canceling"
ECHELON = "echelon"

class SOMStoichiometry(object):

  def __init__(self, som, stoichiometry):
    if not isinstance(som, SOM):
      raise ValueError("First argument must be a SOM.")
    if not util.isFloat(stoichiometry):
      raise ValueError("Second argument must be a float.")
    self.som = som
    self.stoichiometry = stoichiometry
    self.identifier = self.makeId()

  def __repr__(self):
    return self.makeId()

  def makeId(self):
  	return "%s * %2.2f" % (str(self.som), self.stoichiometry)


class SOMReaction(object):

  def __init__(self, reactants, products, label):
    self.reactants = reactants
    self.products = products
    self.label = label
    self.identifier = self.makeId()
    self.category = self.getCategory()

  def __repr__(self):
    return self.makeId()

  def makeId(self):
    """
    Provides a string representation of the som_reaction
    :return str:
    """
    def makeStoichiometryString(som_stoichiometry):
      num = som_stoichiometry.stoichiometry
      if np.isclose(num, 1.0):
        return ''
      else:
        return "%2.2f " % num
    
    def makeTermCollection(som_stoichiometries):
      """
      Formats a set of terms with stoichiometries.
      :param list-MoleculeStoichiometry:
      :return str:
      """
      term_collection = ''
      for s_s in som_stoichiometries:
        term = "%s%s" % (makeStoichiometryString(s_s), str(s_s.som.identifier))
        if len(term_collection) == 0:
          term_collection += term
        else:
          term_collection += " + " + term
      return term_collection
    #
    reactant_collection = makeTermCollection(self.reactants)
    product_collection = makeTermCollection(self.products)
    #
    reaction_str = "%s: %s -> %s" % (self.label,
        reactant_collection, product_collection)
    reaction_str = reaction_str
    return reaction_str

  def getCategory(self):
    """
    Return category of SOMReaction. 
    Return reaction_n_n
    if none of the above applies.
    Each reactant/product (SOMStoichiometry)
    should have a unique SOM wthin the side. 
    :return str:
    """
    num_reactants = len(self.reactants)
    num_products = len(self.products)
    stoichiometry_reactants = [r.stoichiometry for r \
                                  in self.reactants]
    stoichiometry_products = [p.stoichiometry for p \
                             in self.products]
    for reaction_category in cn.REACTION_SUMMARY_CATEGORIES:
      if reaction_category.predicate(num_reactants, num_products, 
                                     stoichiometry_reactants, 
                                     stoichiometry_products):
        return reaction_category.category
    # if none of the above, return reaction_n_n
    return cn.REACTION_n_n


class GAMES_PP(nx.DiGraph):
  """
  Similar to MESGraph, The Message -GAMES++- algorithm creates
  a directed graph of SOMs, and updates the graph using
  the (reduced) row echolon stoichiometry matrix. 
  There are three graphical errors and two matricial errors.

  <Graphical Errors - Type I, II, and III>
  Type I Error occurs when we find inequality between two molecules
  in the same SOM, because each element in a SOM has the same weight.
  Type II Error implies there is cyclism between molecules, such as
  A < B < C < ... < A, which is physically impossible.
  Finally, Type III Errors when it cannot merge two SOMs because
  there is already an existing arc.

  <Matricial Errors - Canceling and Echelon>
  A canceling error occurs when constructing stoichiometry matrix.
  As the algorithm calculates net stoichiometry using SOMs, 
  if a column (a net som_stoichiometry reaction) has only one sign,
  it causes an arror.
  The Echelon error is similar, but it occurs after the matrix is 
  reduced to echelon, or reduced row echelon form. If there is only
  soms with one sign (+ or -), and this causes mass balance error.
  """
  def __init__(self, simple=None):
    """
    :param SimpleSBML simple:
    """
    self.simple = simple
    self.reactions = self._getNonBoundaryReactions(simple)
    self.molecules = self._getNonBoundaryMolecules(simple, self.reactions)
    self.som_stoichiometry_matrix = None
    # reactions before LU decomposition
    self.reactions_lu = []
    # SOMReactinos before LU decomposition
    self.som_reactions_lu = []
    # SOMReactions after LU decomposition
    self.reduced_som_reactions = []
    # RREF SOMReactions after LU -> RREF
    self.rref_som_reactions = []
    # L matrix from LU decomposition
    self.lower = None
    # L^-1 matrix from LU decomposition
    self.lower_inverse  = None
    # P matrix from 'P'LU decomposition
    self.perm_inverse = None
    # permuted stoichiometry matrix
    self.permuted_matrix = None
    # U matrix from LU decomposition
    self.echelon_df = None
    # RREF operation matrix 
    self.rref_operation = None
    # RREF matrix
    self.rref_df = None
    # Components for SOMGraph
    super(GAMES_PP, self).__init__()
    self.soms = self.initializeSOMs(self.molecules)
    # networkx method
    self.add_nodes_from(self.soms)
    self.identifier = self.makeId()
    #
    # List of errors
    # Mass balance error from U matrix (LP decomposition)
    # If not error was detected from U, the error is from RREF
    self.echelon_errors = []
    # Can't add arc
    self.type_one_errors = []
    # Mass balance error from net stoichiometry 
    self.canceling_errors = []
    # SOM cycle ({A} -> {B} -> ... -> {A})
    self.type_two_errors = []
    # Can't merge nodes
    self.type_three_errors = []
    # Cant't add arcs using SOMs (after LU decomposition)
    #self.type_one_som_errors = set()
    # for error return:
    self.error_summary = []
  
  def __repr__(self):
    return self.identifier
  
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
    
    return identifier
  
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

  def convertReactionToSOMReaction(self, reaction):
    """
    Convert simpleSBML Reaction to SOMReaction
    :param Reaction reaction:
    :return SOMReaction:
    """ 
    reactants = reaction.reactants
    products = reaction.products
    reactant_soms = list({self.getNode(r.molecule) for r in reactants})
    product_soms = list({self.getNode(p.molecule) for p in products})
    
    def getSumStoichiometry(som, species):
      """
      :param SOM som:
      :param list-MoleculeStoichiometry species:
      :return SOMStoichiometry:
      """   
      sum_stoichiometry = 0.0
      for s in species:
        if self.getNode(s.molecule) == som:
          sum_stoichiometry += s.stoichiometry
      return SOMStoichiometry(som, sum_stoichiometry)
    #
    ss_reactants = []
    ss_products = []
    for som in reactant_soms:
      ss_reactants.append(getSumStoichiometry(som, reactants))
    for som in product_soms:
      ss_products.append(getSumStoichiometry(som, products))
    #
    return SOMReaction(ss_reactants, ss_products, reaction.label)
  
  def getStoichiometryMatrix(self, reactions, species, som=False):
    """
    Creates a full stoichiometry matrix
    using non-boundary reactions.
    Due to the issue with multiple Molecules with same name, 
    species can also be a list of str - in this case, 
    each component must be molecule.name. 
    Helped by https://gist.github.com/lukauskas/d1e30bdccc5b801d341d
    :param list-Reaction/SOMReaction reactions:
    :param list-str/Molecule/SOM species:
    :return pd.DataFrame:
    """
    reaction_labels = [r.label for r in reactions]
    if som:
      species_names = [s.identifier for s in species]
    elif type(species[0]) == str:
      species_names = species
    else:
      species_names = [m.name for m in species]
    stoichiometry_matrix = pd.DataFrame(
        0.0,
        index=species_names,
        columns=reaction_labels
        )
    for reaction in reactions:
      if som:
        reactants = {r.som.identifier:r.stoichiometry for r in reaction.reactants}
        products = {p.som.identifier:p.stoichiometry for p in reaction.products}
      else:
        reactants = {r.molecule.name:r.stoichiometry for r in reaction.reactants}
        products = {p.molecule.name:p.stoichiometry for p in reaction.products}
      reaction_species = list(set(reactants.keys()).union(products.keys()))
      for species_name in reaction_species:
        net_stoichiometry = products.get(species_name, 0.0) - reactants.get(species_name, 0.0)
        stoichiometry_matrix.loc[species_name, reaction.label] = net_stoichiometry
    return stoichiometry_matrix
  
  def decomposeMatrix(self, mat_df):
    """
    LU decomposition of the stoichiometry matrix.
    First it transposes the input matrix
    and find L, U matrices and permutation list. 
    :param pandas.DataFrame mat_df:
    :return pandas.DataFrame echelon_df:
    """
    mat_t = mat_df.T
    idx_mat_t = mat_t.index
    cols_mat_t = mat_t.columns
    #
    diff = None
    if mat_t.shape[0] > mat_t.shape[1]:
      diff = mat_t.shape[0] - mat_t.shape[1]
      for i in range(diff):
        mat_t["_" + str(i)] = np.zeros(mat_t.shape[0])
    # LU decomposition
    perm, lower, upper = lu(mat_t)
    # #### Trying to round up lower and upper, to avoid precision issue related to 0.0
    # lower = np.round(lower_raw, 3)
    # upper = np.round(upper_raw, 3)
    perm_inverse = perm.T
    permuted_m = (perm_inverse).dot(mat_t)
    pivot_index = [list(k).index(1) for k in perm_inverse]
    new_idx_mat_t = [idx_mat_t[idx] for idx in pivot_index]
    # we save as; lower_inverse * perm_df = echelon_df
    # if we added zero columns previously, delete them. 
    if diff:
      perm_df = pd.DataFrame(permuted_m,
          index=new_idx_mat_t,
          columns=mat_t.columns).drop(columns = mat_t.columns[-diff:]).T
      echelon_df = pd.DataFrame(upper,
          index=new_idx_mat_t,
          columns=mat_t.columns).drop(columns = mat_t.columns[-diff:]).T
    else:
      perm_df = pd.DataFrame(permuted_m,
          index=new_idx_mat_t,
          columns=mat_t.columns).T
      echelon_df = pd.DataFrame(upper,
          index=new_idx_mat_t,
          columns=mat_t.columns).T
    lower_inverse = pd.DataFrame(inv(lower),
          index=new_idx_mat_t,
          columns=new_idx_mat_t)
    self.perm_inverse = perm_inverse
    self.permuted_matrix = perm_df
    self.lower = lower
    self.lower_inverse = lower_inverse
    self.echelon_df = echelon_df
    return echelon_df
  
  def getRREFMatrix(self, echelon_df):
    """
    Get RREF of the stoichiometry matrix.
    Stoichiometry matrix should be in echelon form.
    Columns are reactions and indices are species. 
    :param pandas.DataFrame echelon_df:
    :return pandas.DataFrame rref_df:
    """
    rref_operation = np.identity(echelon_df.T.shape[0])
    rref_operation = pd.DataFrame(rref_operation,
                     index = echelon_df.columns,
                     columns = echelon_df.columns)

    # now update the operation matrix, finally multiply (dot product) two matrices
    for idx, colname in enumerate(echelon_df.columns):
      reaction_series = echelon_df[colname]
      # Find the first nonzero values
      # Deprecation: instead of np.nonzero(Series), use Series.to_numpy().nonzero()
      # to fix error on GitHub - not use to_numpy here
      ##nonzero_idx = echelon_df[colname].to_numpy().nonzero()[0]
      nonzero_idx = np.array([idx for idx, val in enumerate(echelon_df[colname]) if val != 0])
      # Skip if there is no nonzero value or if it is first reaction
      if not nonzero_idx.any() or idx == 0:
        continue
      nonzero_species = reaction_series.index[nonzero_idx[0]]
      nonzero_value = reaction_series[nonzero_idx[0]]
      # find current nonzero index
      current_echelon_t = rref_operation.dot(echelon_df.T).T
      for prev_colname in current_echelon_t.columns[:idx]:
        if np.round(current_echelon_t[prev_colname][nonzero_species], 3) != 0.0:
          reduction_value = current_echelon_t[prev_colname][nonzero_species]
          rref_operation.at[prev_colname, colname] = (-1.0) * reduction_value / nonzero_value
    rref_df = np.round(rref_operation.dot(echelon_df.T).T, 3)
    self.rref_operation = rref_operation	
    self.rref_df = rref_df
    return rref_df

  def convertMatrixToSOMReactions(self, mat_df):
    """
    Convert a stoichiometry matrix to SOMReactions,
    where columns are SOMReactions and 
    rows are SOMs (species).
    :param pandas.DataFrame mat_df:
    :return list-SOMReaction reactions:
    """
    reactions = []
    for reaction_name in mat_df.columns:
      reaction_elements = mat_df[reaction_name]
      reactants = [SOMStoichiometry(
          self.getNode(som_label),
          np.round(abs(reaction_elements[som_label]), 3)
          ) \
          for som_label in reaction_elements.index \
          if reaction_elements[som_label]<TOLERANCE*(-1)]
      products = [SOMStoichiometry(
          self.getNode(som_label),
          np.round(abs(reaction_elements[som_label]), 3)
          ) \
          for som_label in reaction_elements.index \
          if reaction_elements[som_label]>TOLERANCE]
      reactions.append(SOMReaction(
          reactants=reactants,
          products=products,
          label=reaction_name
          ))
    return reactions
  
  def getNode(self, input_arg):
    """
    Find a node(SOM) containing the given molecule
    or the SOM identifier. 
    If no such SOM exists, return False
    :param Molecule/str input_arg:
    :return SOM/False:
    """
    for som in list(self.nodes):
      if isinstance(input_arg, Molecule):
        for mole in som.molecules:
          if mole.name == input_arg.name:
            return som
      elif isinstance(input_arg, str):
      	if som.identifier == input_arg:
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
  
  def addReaction(self, reaction=None):
    """
    Add a reaction to self.reactions_lu
    :param reaction Reaction:
    """
    if reaction not in self.reactions_lu:
      self.reactions_lu.append(reaction)  
  
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
  
  def addArc(self, arc_source, arc_destination, reaction):
    """
    Add a single arc (edge) using two SOMs and reaction/somreaction.
    :param SOM arc_source:
    :param SOM arc_destination:
    :param Reaction/SOMReaction reaction:
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
      error_count = 0
      for arc in arcs:
        if not self.checkTypeOneError(arc, reaction):
          som_source = self.getNode(arc[0])
          som_destination = self.getNode(arc[1])
          self.addArc(som_source, som_destination, reaction)
        else:
          error_count += 1
      # the reaction has not added any errors, will be used for lu
      if error_count == 0:
        self.reactions_lu.append(reaction)
      self.identifier = self.makeId()
  
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
      error_count = 0
      for arc in arcs:
        if not self.checkTypeOneError(arc, reaction):
          som_source = self.getNode(arc[0])
          som_destination = self.getNode(arc[1])
          self.addArc(som_source, som_destination, reaction)
        else:
          error_count += 1
      # the reaction has not added any errors, will be used for lu
      if error_count == 0:
        self.reactions_lu.append(reaction)
      self.identifier = self.makeId()
  
  def processEqualSOMReaction(self, reaction):
    """
    Process a 1-1 SOMReaction to check 
    mergeability of nodes.
    If there are existing arcs, two nodes
    cannot be merged which will lead to a 
    type_three_error.
    :param SOMReaction reaction:
    """
    if reaction.category != cn.REACTION_1_1:
      pass
    else:
      reactant = reaction.reactants[0].som
      product = reaction.products[0].som
      if self.has_edge(reactant, product) or self.has_edge(product, reactant):
        self.type_three_errors.append(reaction)

  def processUnequalSOMReaction(self, reaction):
    """
    Process a 1-n or n-1 SOMReaction to add arcs.
    An arc flows from SOM with less mass (source)
    to SOM with greater mass (destination). 
    :param SOMReaction reaction:
    """
    category = reaction.category
    if (category!=cn.REACTION_n_1) and (category!=cn.REACTION_1_n):
      pass
    else:
      if category == cn.REACTION_n_1:
        destination = [reaction.products[0].som]
        source = [reactant.som for reactant in reaction.reactants]
      elif category == cn.REACTION_1_n:
        destination = [reaction.reactants[0].som]
        source = [product.som for product in reaction.products]
      arcs = itertools.product(source, destination)
      for arc in arcs:
        if arc[0] == arc[1]:
          #self.type_one_som_errors = self.type_one_som_errors.add(reaction.label)
          raise ValueError("SOM Unequal reaction leads an arc between two soms")
        else:
          self.addArc(arc[0], arc[1], reaction)
      self.identifier = self.makeId()
  
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
  
  def checkTypeOneError(self, arc, inequality_reaction):
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
  
  # We may remove this method
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
        # We are not using addTypeTwoError yet. 
        # self.addTypeTwoError(cycle)
        self.type_two_errors.append(cycle)
        # if not self.type_two_error:
        #   self.type_two_error = True
      return True
  
  def processErrorReaction(self, reaction):
    """
    Simply add error reactions to error
    :param Reaction reaction:
    :return bool:
    """
    error_reactions = [r.label for r in self.echelon_errors]
    if reaction.label not in error_reactions:
      self.echelon_errors.append(reaction)
    ## Here, we may need to add more info that will 
    ## help us track the operations that lead to this error 
    return True
  
  def analyze(self, reactions=None, simple_games=False, rref=True, error_details=False, suppress_message=False):
    """
    Using the stoichiometry matrix, compute
    row reduced echelon form and create SOMGraph
    Add arcs or sending error messages using
    checkTypeOneError or checkTypeTwoError.
    :param list-Reaction reactions:
    :param bool rref:
    :param bool error_details:
    :param bool suppress_message:
    :return bool:
    """
    multimulti_error_found = False
    if reactions is None:
      reactions = self.simple.reactions
    # Associate the reaction category with the function
    # that processes that category
    report = cn.NULL_STR
    reaction_dic = {
        cn.REACTION_1_1: self.processUniUniReaction,
        cn.REACTION_1_n: self.processUniMultiReaction,
        cn.REACTION_n_1: self.processMultiUniReaction,
        cn.REACTION_n_n: self.addReaction,
        }
    # Process each type of reaction - Type I error will be detected here
    for category in reaction_dic.keys():
      for reaction in [r for r in reactions if r.category == category]:
        func = reaction_dic[category]
        func(reaction)
    # detect type II error
    self.checkTypeTwoError()
    #########################
    # if we find type I or II errors, we make it a simple_game
    if self.type_one_errors or self.type_two_errors:
      simple_games = True
    ########################
    # if simple_games, we only run elementary operations
    if not simple_games:
      # reaction_lu are reactions for LU decomposition, i.e. multi-multi reactions
      if self.reactions_lu:
        for reaction in self.reactions_lu:
          self.som_reactions_lu.append(
              self.convertReactionToSOMReaction(reaction)
              )
        # Now, step 1: creates SOMStoichiometryMatrix
        self.som_stoichiometry_matrix = self.getStoichiometryMatrix(self.som_reactions_lu, list(self.nodes), som=True)
        # step 2: reconvert it into SOMReactions and examine 'canceling errors'
        initial_reduced_som_reactions = self.convertMatrixToSOMReactions(self.som_stoichiometry_matrix)
        dropping_reactions = []
        for r in initial_reduced_som_reactions:
          if r.category == cn.REACTION_ERROR:
            self.canceling_errors.append(r)
            dropping_reactions.append(r.label)
            multimulti_error_found = True
        # step 3: decompose using LU decompositon and check errors (echelon, type_three)
        if not multimulti_error_found:
          echelon_df = self.decomposeMatrix(self.som_stoichiometry_matrix)
          self.reduced_som_reactions = self.convertMatrixToSOMReactions(echelon_df)
          som_reaction_dic = {
              cn.REACTION_ERROR: self.processErrorReaction,
              cn.REACTION_1_1: self.processEqualSOMReaction,
              cn.REACTION_1_n: self.processUnequalSOMReaction,
              cn.REACTION_n_1: self.processUnequalSOMReaction,
              }
          for category in som_reaction_dic.keys():
            for reaction in [r for r in self.reduced_som_reactions if r.category == category]:
              func = som_reaction_dic[category]
              func(reaction)
          # checking if there was any error by LU decompoistion
          if self.echelon_errors or self.type_three_errors:
            multimulti_error_found = True
        # step 4: get RREF and check errors (same as LU decomposition case)
        if not multimulti_error_found:
          rref_df = self.getRREFMatrix(self.echelon_df)
          self.rref_som_reactions = self.convertMatrixToSOMReactions(rref_df)
          for category in som_reaction_dic.keys():
            for reaction in [r for r in self.rref_som_reactions if r.category == category]:
              func = som_reaction_dic[category]
              func(reaction)
    if not suppress_message:
      print("Model analyzed...")
    if error_details:
      print("Type I error: ", self.type_one_errors)
      print("Type II error: " , self.type_two_errors)
      print("Canceling error: ", self.canceling_errors)
      print("Echelon error: ", self.echelon_errors)
      print("Type III error: ", self.type_three_errors, "\n")
    if self.echelon_errors or self.type_one_errors or self.type_two_errors \
        or self.canceling_errors or self.type_three_errors:
      if self.type_one_errors:
        self.error_summary.append(ErrorSummary(type=TYPE_I, errors=self.type_one_errors))
      if self.type_two_errors:
        self.error_summary.append(ErrorSummary(type=TYPE_II, errors=self.type_two_errors))
      if self.type_three_errors:
        self.error_summary.append(ErrorSummary(type=TYPE_III, errors=self.type_three_errors))
      if self.canceling_errors:
        self.error_summary.append(ErrorSummary(type=CANCELING, errors=self.canceling_errors))
      if self.echelon_errors:
        self.error_summary.append(ErrorSummary(type=ECHELON, errors=self.echelon_errors))
      if not suppress_message:
        print("At least one error found.\n")
      return True
    else:
      if not suppress_message:
        print("No error found.")
      return False







