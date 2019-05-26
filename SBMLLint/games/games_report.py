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
    Generate report and store it in 
    self.report_type_one_errors.
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
      report = report + "*"*40 + "\n"
    report = report + "-"*50 + "\n"
    return report

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
  	return linked_molecules, linked_reactions

  def reportEchelonError(self, echelon_errors):
    """
    Generate report for echelon errors, i.e.
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
    linkage_report = NULL_STR
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
      print("reported_reactions", reported_reactions)
      # part 2: find reactions that created related SOMs
      nonzero_elements = self.mesgraph.som_stoichiometry_matrix.index
      for r in reported_reactions:
        som_row = self.mesgraph.som_stoichiometry_matrix[r]
        nonzero_elements = nonzero_elements.intersection(
        	som_row[som_row.to_numpy().nonzero()[0]].index
        	)
      canceled_soms = nonzero_elements.difference(nonzero_result_series.index)
      print("canceled_soms", canceled_soms)
      for som_name in canceled_soms:
        print("som_name is ", som_name)
        linked_molecules, linked_reactions = self.getMoleculeLinkage(som_name, reported_reactions)
        linkage_report = linkage_report + "The following molecules,\n"
        for m in list(linked_molecules):
          linkage_report = linkage_report + m + ", "
        linkage_report = linkage_report + "\nHave equal mass by the following reactions.\n"
        for r in linked_reactions:
          reaction = self.mesgraph.simple.getReaction(r)
          linkage_report = linkage_report + reaction.makeIdentifier(is_include_kinetics=False)
          linkage_report = linkage_report + "\n"
      #
      # generate an error report for a single echelon error
      print("linkage report :")
      print(linkage_report)
      error_report = "The following reactions create a mass imbalance.\n\n"
      for r in reported_reactions:
      	simple_reaction = self.mesgraph.simple.getReaction(r)
      	error_report = error_report + simple_reaction.makeIdentifier(is_include_kinetics=False)
      	error_report = error_report + "\n"

      echelon_report = echelon_report + "\n" + error_report + "*"*40 + "\n"
    echelon_report = echelon_report + "-"*50 + "\n"
    return echelon_report








