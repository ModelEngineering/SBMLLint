"""
Reporting class for GAMES Plus (GAMES_PP) algorithm
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common.reaction import Reaction
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.games.games_pp import SOMStoichiometry, SOMReaction, GAMES_PP
from SBMLLint.games.som import SOM
from SBMLLint.common import simple_sbml

import networkx as nx
import numpy as np
import pandas as pd
import os
import re
import tesbml

NULL_STR = ""


class GAMESReport(object):

  def __init__(self, mesgraph, errors=None):
  	self.mesgraph = mesgraph
    self.errors = errors
    self.report_type_one_errors = NULL_STR

  def getSOMPath(self, som, mole1, mole2):
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
        # if reaction.label is not already included in the attribute,
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

  def buildSOMPathReport(self, molecule_name1, molecule_name2):
    """
    Print out shortest path between two molecules within a 
    same SOM. Molecule names (str) are passed as arguments. 
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
      # in case molecule_name1 == molecule_name2
      if molecule_name1 == molecule_name2:
        path_report = path_report + "Clearly, %s %s %s\n" % (
            molecule_name1, cn.EQUAL, molecule_name2)
      else:
        som_path = self.mesgraph.getSOMPath(som1, 
                                   self.mesgraph.simple.getMolecule(molecule_name1), 
                                   self.mesgraph.simple.getMolecule(molecule_name2))
        for pat in som_path:
          path_report = path_report + "\n%s %s %s by reaction(s):\n" % (pat.node1, cn.EQUAL, pat.node2)
          for r in pat.reactions:
            som_reaction = self.mesgraph.simple.getReaction(r)
            path_report = path_report + "%s\n" % (som_reaction.makeIdentifier(is_include_kinetics=False))
      return path_report

  def reportTypeOneError(self):
    """
    Generate report and store it in 
    self.report_type_one_errors.
    :return: str type_one_report
    """
  	type_one_report = NULL_STR
  	pass

