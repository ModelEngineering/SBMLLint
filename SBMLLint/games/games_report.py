"""
Reporting class for GAMES Plus (GAMES_PP) algorithm
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common.reaction import Reaction
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.games.games_pp import SOMStoichiometry, SOMReaction, GAMES_PP, TOLERANCE
from SBMLLint.games.som import SOM
from SBMLLint.common import simple_sbml

import collections
import networkx as nx
import numpy as np
import pandas as pd
import tesbml

NULL_STR = ""
ReactionOperation = collections.namedtuple("ReactionOperation", 
    "reaction operation")
NUM_STAR = 50


class GAMESReport(object):

  def __init__(self, mesgraph, errors=None):
  	self.mesgraph = mesgraph
  	self.errors = errors
  	self.report_type_one_errors = NULL_STR

  def getMoleculeEqualityPath(self, som, mole1, mole2):
    """
    Create an undirected graph between
    two molecules within a SOM
    and find the shortest path
    :param SOM som:
    :param Molecule mole1:
    :param Molecule mole2:
    :return PathComponents: som_path
    """   
    molecule1 = mole1.name
    molecule2 = mole2.name
    # construct undirected graph
    subg = nx.Graph()
    # here, every reaction is 1-1 reaction
    for reaction in list(som.reactions):
      node1 = reaction.reactants[0].molecule.name
      node2 = reaction.products[0].molecule.name
      if subg.has_edge(node1, node2):
        reaction_label = subg.get_edge_data(node1, node2)[cn.REACTION]
        # if reaction.label is not included in the attribute, add its label
        if reaction.label not in set(reaction_label):
          reaction_label = reaction_label + [reaction.label]
      else:
        reaction_label = [reaction.label]    
      subg.add_edge(node1, node2, reaction=reaction_label)
    path = [short_p for short_p in nx.shortest_path(subg, 
                                                    source=molecule1, 
                                                    target=molecule2
                                                    )
           ]
    som_path = []
    for idx in range(len(path)-1):
      edge_reactions = subg.get_edge_data(path[idx], path[idx+1])[cn.REACTION]
      som_path.append(cn.PathComponents(node1=path[idx], 
                                        node2=path[idx+1],
                                        reactions=edge_reactions
                                        )
                     )
    return som_path

  def getMoleculeEqualityPathReport(self, molecule_name1, molecule_name2):
    """
    Print out shortest path between two molecules within a 
    same SOM. Molecule names (str) are arguments. 
    :param str molecule_name1:
    :param str molecule_name2:
    :return bool/str: path_report
    """
    path_report = NULL_STR
    som1 = self.mesgraph.getNode(self.mesgraph.simple.getMolecule(molecule_name1))
    som2 = self.mesgraph.getNode(self.mesgraph.simple.getMolecule(molecule_name2))
    if som1 != som2:
      return False
    else:
      # in case molecule_name1 == molecule_name2, for example A -> A + B
      if molecule_name1 == molecule_name2:
        path_report = path_report + "Clearly, %s %s %s\n" % (
            molecule_name1, cn.EQUAL, molecule_name2)
      else:
        som_path = self.getMoleculeEqualityPath(som1, 
                                   self.mesgraph.simple.getMolecule(molecule_name1), 
                                   self.mesgraph.simple.getMolecule(molecule_name2))
        for pat in som_path:
          path_report = path_report + "\n%s %s %s by reaction(s):\n" % (pat.node1, cn.EQUAL, pat.node2)
          for r in pat.reactions:
            som_reaction = self.mesgraph.simple.getReaction(r)
            path_report = path_report + "%s\n" % (som_reaction.makeIdentifier(is_include_kinetics=False))
      return path_report

  def getMoleculeInequalityPathReport(self, molecule_name1, molecule_name2, reaction_names):
  	"""
  	Print the reactions that infer inequality between molecules. 
  	Molecule name are given as arguments.
  	Especially, the mass of molecule 1 is less than that of molecule 2
  	:param str molecule_name1:
  	:param str molecule_name2:
  	:param list-str reaction_names:
  	:return str: path_report
  	"""
  	path_report = NULL_STR
  	path_report = path_report + "%s %s %s by reaction(s):\n" % (molecule_name1, cn.LESSTHAN, molecule_name2)
  	for reaction_name in reaction_names:
  	  reaction = self.mesgraph.simple.getReaction(reaction_name)
  	  path_report = path_report + "%s\n" % (reaction.makeIdentifier(is_include_kinetics=False))
  	return path_report

  def reportTypeOneError(self, type_one_errors):
    """
    Generate report for Type I Errors. 
    Type I error occurs when there is a reaction
    that needs to add an arc between molecules,
    while the molecules are already included
    in the same SOM.
    :param list-PathComponents type_one_errors:
    :return str: type_one_report
    """
    report = NULL_STR
    for pc in type_one_errors:
      mole1 = pc.node1
      mole2 = pc.node2
      reactions = pc.reactions
      report = report + self.getMoleculeEqualityPathReport(mole1, mole2)
      report = report + "\nHowever, "
      report = report + self.getMoleculeInequalityPathReport(mole1, mole2, reactions)
      report = report + "*"*NUM_STAR + "\n"
    report = report + "-"*NUM_STAR + "\n"
    return report

  def decompostSOMCycle(self, cycle):
  	"""
  	Create list of PathComponents
  	using a cycle of SOMs
  	:param list-SOM cycle:
  	:return list-PathComponents: error_cycle
  	"""
  	error_cycle = []
  	cycle2 = cycle[1:] + [cycle[0]]
  	for first, second in zip(cycle, cycle2):
  	  som1_moles = {m.name for m in list(first.molecules)}
  	  som2_moles = {m.name for m in list(second.molecules)}
  	  reaction_data = self.mesgraph.get_edge_data(first, second)[cn.REACTION]
  	  nodes1 = []
  	  nodes2 = []
  	  reaction_labels = []  
  	  for r in reaction_data:
  	  	reaction = self.mesgraph.simple.getReaction(r)
  	  	if reaction.category == cn.REACTION_n_1:
  	  	  sources = {r.molecule.name for r in reaction.reactants}
  	  	  destinations = {p.molecule.name for p in reaction.products}
  	  	elif reaction.category == cn.REACTION_1_n:
  	  	  sources = {p.molecule.name for p in reaction.products}
  	  	  destinations = {r.molecule.name for r in reaction.reactants}
  	  	node2 = list(destinations.intersection(som2_moles))[0]
  	  	for node1 in list(sources.intersection(som1_moles)):
  	  	  nodes1.append(node1)
  	  	  nodes2.append(node2)
  	  	  reaction_labels.append(reaction.label)
  	  error_cycle.append(cn.PathComponents(node1=nodes1,
  	  	                                   node2=nodes2,
  	  	                                   reactions=reaction_labels))
  	return error_cycle

  def getSOMPath(self, som, molecule1, molecule2):
  	"""
  	Create an undirected graph 
  	between two molecules within a SOM
  	and find the shortest path
  	:param SOM som:
  	:param str molecule1:
  	:param str molecule2:
  	:return list-PathComponents: som_path
  	"""
  	# undirected graph
  	subg = nx.Graph()
  	for reaction in list(som.reactions):
  	  node1 = reaction.reactants[0].molecule.name
  	  node2 = reaction.products[0].molecule.name
  	  if subg.has_edge(node1, node2):
  	    reaction_label = subg.get_edge_data(node1, node2)[cn.REACTION]
  	    # if reaction.label is not already included in the attribute
  	    if reaction.label not in set(reaction_label):
  	      reaction_label = reaction_label + [reaction.label]
  	  else:
  	    reaction_label = [reaction.label]
  	  subg.add_edge(node1, node2, reaction=reaction_label)
  	path = [short_p for short_p in nx.shortest_path(subg,
  	                                                source=molecule1,
  	                                                target=molecule2)]
  	som_path = []
  	for idx in range(len(path)-1):
  	  edge_reactions = subg.get_edge_data(path[idx], path[idx+1])[cn.REACTION]
  	  som_path.append(cn.PathComponents(node1=path[idx],
  	                                    node2=path[idx+1],
  	                                    reactions=edge_reactions))
  	return som_path

  def reportTypeTwoError(self, type_two_errors):
  	"""
  	Generate report for Type II Errors.
  	Type II Error occurs when there is
  	a cycle between SOMs, which
  	should not happen.
  	:param list-(list-SOMs) type_two_errors:
  	:return str: type_two_report
  	"""
  	report = NULL_STR
  	for cycle in type_two_errors:
  	  error_cycle = self.decompostSOMCycle(cycle)
  	  last_node = None
  	  for one_path in error_cycle:
  	  	comb = zip(one_path.node1, one_path.node2, one_path.reactions)
  	  	som_path_report = NULL_STR
  	  	if last_node is None:
  	  	  first_molecule = one_path.node1[0]
  	  	  last_node = one_path.node2[0]
  	  	else:
  	  	  if last_node in one_path.node1:
  	  	  	pass
  	  	  else:
  	  	    som = self.mesgraph.getNode(self.mesgraph.simple.getMolecule(last_node))
  	  	    som_path = self.getSOMPath(som, last_node, one_path.node1[0])
  	  	    som_path_report = som_path_report + self.reportSOMPath(som_path)

  	  	for tup in comb:
  	  	  report = report + "\n%s %s %s by\n" % (tup[0], cn.LESSTHAN, tup[1])
  	  	  reaction = self.mesgraph.simple.getReaction(tup[2])
  	  	  report = report + "%s\n" % (reaction.makeIdentifier(is_include_kinetics=False))
  	  	if len(one_path.node1) > 1:
  	  	  for node1, node2 in zip(one_path.node1, one_path.node1[1:]):
  	  	    mole1 = self.mesgraph.simple.getMolecule(node1)
  	  	    mole2 = self.mesgraph.simple.getMolecule(node2)
  	  	    som = self.mesgraph.getNode(mole1)
  	  	    som_path = self.getSOMPath(som, node1, node2)
  	  	    report = report + self.reportSOMPath(som_path)
  	  	report = report + som_path_report
  	  last_molecule = tup[1]
  	  if first_molecule != last_molecule:
  	  	som = self.mesgraph.getNode(self.mesgraph.simple.getMolecule(first_molecule))
  	  	som_path = self.getSOMPath(som, first_molecule, last_molecule)
  	  	report = report + self.reportSOMPath(som_path)
  	  report = report + "*"*NUM_STAR + "\n\n"
  	report = report + "\n" + "-"*NUM_STAR + "\n"
  	return report

  def reportSOMPath(self, som_path):
  	"""
  	Generate a path report between molecules.
  	:param list-PathComponents som_path:
  	:return str: path_report
  	"""
  	path_report = NULL_STR
  	# path_report = path_report + "\nThe equality between molecules are as follows:\n"
  	for one_path in som_path:
  	  path_report = path_report + "\n%s %s %s by\n" % (one_path.node1, cn.EQUAL, one_path.node2)
  	  for r in one_path.reactions:
  	  	reaction = self.mesgraph.simple.getReaction(r)
  	  	path_report = path_report + "%s\n" % (reaction.makeIdentifier(is_include_kinetics=False))
  	return path_report

  def convertOperationSeriesToReactionOperations(self, operation):
    """
    Convert an operation series to 
    a list of reaction operations.
    An 'operation' is either a row or column 
    in an operation matrix, 
    where both column and row indices are reactions. 
    :param pandas.Series operation:
    :return list-ReactionOperation: operations
    """
    operations = []
    nonzero_idx = operation.to_numpy().nonzero()[0]
    nonzero_op = operation[nonzero_idx]
    for idx in range(len(nonzero_op)):
      reaction_op = ReactionOperation(reaction=nonzero_op.index[idx],
      	                              operation=nonzero_op[idx]
      	                              )
      operations.append(reaction_op)
    return operations

  def findSOM(self, som_name):
  	"""
  	Return a SOM by a SOM identifier.
  	:param str som_name:
  	:return SOM/None:
  	"""
  	for node in self.mesgraph.nodes:
  	  if node.identifier == som_name:
  	    return node
  	return None

  def getMoleculeLinkage(self, som_name, reactions):
  	"""
  	Create two lists. 
  	1. molecules in the reported_reactions that are in the same som
  	2. reactions used to merge the molecules
  	:param str som_name:
  	:param list-str reactions:
  	:return list-str: linked_molecules
  	:return list-str: linked_reactions
  	"""
  	pass
  	som = self.findSOM(som_name)
  	molecules = {m.name for m in som.molecules}
  	# linked_molecules: molecules within both the SOM and given reactions
  	linked_molecules = set()
  	for r in reactions:
  	  reaction = self.mesgraph.simple.getReaction(r)
  	  reactants = {m.molecule.name for m in reaction.reactants}
  	  products = {m.molecule.name for m in reaction.products}
  	  som_molecules = molecules.intersection(reactants)
  	  som_molecules = som_molecules.union(molecules.intersection(products))
  	  linked_molecules = linked_molecules.union(som_molecules)
  	linked_reactions = []
  	for sr in som.reactions:
  	  sreactants = {m.molecule.name for m in sr.reactants}
  	  sproducts = {m.molecule.name for m in sr.products}
  	  if sreactants.intersection(linked_molecules) and \
  	      sproducts.intersection(linked_molecules):
  	    linked_reactions.append(sr.label)
  	return list(linked_molecules), linked_reactions

  def reportLinkage(self, linked_molecules, linked_reactions):
  	"""
  	Generate a report for linked molecules.
  	Molecules are linked within a SOM
  	by appropriate reations. 
  	:param list-str linked_molecules:
  	:param list-str linked_reactions:
  	:return str: linkage_report
  	"""
  	linkage_report = NULL_STR
  	linkage_report = linkage_report + "\n" + "->"*int((NUM_STAR/2))
  	if len(linked_molecules) == 1:
  	  m = linked_molecules[0]
  	  linkage_report = linkage_report + "\n%s is a common element in the reactions above.\n" % m
  	else:
  	  linkage_report = linkage_report + "\nThe following molecules,\n"
  	  for m in linked_molecules:
  	    linkage_report = linkage_report + m + "\n"
  	  linkage_report = linkage_report + "Have equal mass by the following reaction(s).\n"
  	  for r in linked_reactions:
  	  	reaction = self.mesgraph.simple.getReaction(r)
  	  	linkage_report = linkage_report + reaction.makeIdentifier(is_include_kinetics=False)
  	  	linkage_report = linkage_report + "\n"
  	linkage_report = linkage_report + "<-"*int((NUM_STAR/2)) + "\n"
  	return linkage_report

  def reportEchelonError(self, echelon_errors):
    """
    Generate a report for echelon errors, i.e.
    mass balance errors from LU decomposition/RREF
    and store in self.report_echelon_errors.
    The operation_df is either mesgraph.lower_inverse
    if LU decomposition yielded echelon_errors,
    and rref_df.dot(lower_inverse) if RREF created errors. 
    :param list-SOMReaction echelon_errors:
    :return str: echelon_report
    """
    if self.mesgraph.rref_df is None:
      operation_df = self.mesgraph.lower_inverse
    else:
      operation_df = self.mesgraph.rref_operation.dot(self.mesgraph.lower_inverse)
    echelon_report = NULL_STR
    error_report = NULL_STR
    for reaction in echelon_errors:
      reaction_label = reaction.label
      operation_series = operation_df.T[reaction_label]
      if self.mesgraph.rref_df is None:
        result_series = self.mesgraph.echelon_df[reaction_label]
      else:
        result_series = self.mesgraph.rref_df[reaction_label]
      nonzero_result_series = result_series[result_series.to_numpy().nonzero()[0]]
      # part 1: find reactions that caused mass balance errors
      reaction_operations = self.convertOperationSeriesToReactionOperations(operation_series)
      reported_reactions = [r.reaction for r in reaction_operations]
      # part 2: find reactions that created related SOMs
      nonzero_elements = self.mesgraph.som_stoichiometry_matrix.index
      for r in reported_reactions:
        som_row = self.mesgraph.som_stoichiometry_matrix[r]
        nonzero_elements = nonzero_elements.intersection(
        	som_row[som_row.to_numpy().nonzero()[0]].index
        	)
      canceled_soms = nonzero_elements.difference(nonzero_result_series.index)
      linkage_report = NULL_STR
      for som_name in canceled_soms:
        linked_molecules, linked_reactions = self.getMoleculeLinkage(som_name, reported_reactions)
        linkage_report = linkage_report + self.reportLinkage(linked_molecules, linked_reactions)
        # linkage_report = linkage_report + "\n" + "->"*int((NUM_STAR/2))
        # if len(linked_molecules) == 1:
        #   m = list(linked_molecules)[0]
        #   linkage_report = linkage_report + "\n%s is a common element in the reactions above.\n" % m
        # else:
        #   linkage_report = linkage_report + "\nThe following molecules,\n"
        #   for m in list(linked_molecules):
        #     linkage_report = linkage_report + m + "\n"
        #   linkage_report = linkage_report + "Have equal mass by the following reaction(s).\n"
        #   for r in linked_reactions:
        #     reaction = self.mesgraph.simple.getReaction(r)
        #     linkage_report = linkage_report + reaction.makeIdentifier(is_include_kinetics=False)
        #     linkage_report = linkage_report + "\n"
        # linkage_report = linkage_report + "<-"*int((NUM_STAR/2)) + "\n"
      # generate an error report for a single echelon error
      error_report = "The following reactions create a mass imbalance:\n\n"
      for r in reported_reactions:
      	simple_reaction = self.mesgraph.simple.getReaction(r)
      	error_report = error_report + simple_reaction.makeIdentifier(is_include_kinetics=False)
      	error_report = error_report + "\n"
      echelon_report = echelon_report + "\n" + error_report + linkage_report + "\n" + "*"*NUM_STAR + "\n"
    echelon_report = echelon_report + "-"*NUM_STAR + "\n"
    return echelon_report








