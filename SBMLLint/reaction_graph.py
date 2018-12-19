# -*- coding: utf-8 -*-
"""Doctrings on reaction_graph.py.

The ReactionGraph class represents a directed graph class between reactants and products
in each reaction. Each reactant is connected to a reaction node (identified by reaction name)
which is then connected to each product. 

Example:
    You should already have an SBML model loaded, which we call 'model'.

    >> import reaction_graph as rg
    >> rg_class = rg.ReactionGraph(model)
    >> rg_class.buildReactionGraph()

"""


import constants as cn
import stoichiometry_matrix as sm

import itertools  
import networkx as nx
import os
import tellurium as te
import tesbml


class ReactionGraph():
    """Creates a reaction graph. 

    This class takes one input argument, an SBML model, and creates
    a directed graph (reaction_graph). The user can use a method to check consistency. 

    Attributes:
        model (SBML model): An SBML model. The argument needed to create a ModelGraph object. 
        consistent (bool): True if model is consistent, False if not. 
        reaction_graph (networkx.DiGraph) : A directed graph between reactants-reaction_name-products
        inconsistent_set (set) : A set of reactions that are inconsistent; the remaining reactions 
                                 within the model should be consistent. 


    """    
    def __init__(self, model):
        self.model = model
        self.consistent = None
        self.reaction_graph = None
        self.inconsistent_set = None

        # if not isinstance(self.model, tesbml.libsedml.Model):
        if self.model == None:
            raise TypeError("Model doesn't exist")

    def buildReactionGraph(self):
        """Creates a reaction graph from a model.
        
        We assume each reaction has at least one reactant and one product.
        if either is undefined as EmptySet, we create a new boundary species 
        whose name is created by reaction name + _BDRY_RCT(_PDT if product). 
       
        Args:
            None. 

        Returns:
            networkx.classes.digraph.DiGraph    
        """
        reaction_list = [reaction.getId() for reaction in self.model.getListOfReactions()]

        G = nx.DiGraph()
        for reaction_name in reaction_list:
            reaction = self.model.getReaction(reaction_name)
            reactant_list = [a.getSpecies() for a in reaction.getListOfReactants()]
            if (reactant_list == []) | (reactant_list == [cn.EMPTYSET]):
                reactant_list = [reaction_name + cn.BDRY_RCT]
            product_list = [a.getSpecies() for a in reaction.getListOfProducts()]
            if (product_list == []) | (product_list == [cn.EMPTYSET]):
                product_list = [reaction_name + cn.BDRY_PDT]
            
            # Add edges to the DiGraph object
            G.add_edges_from(itertools.product(reactant_list, [reaction_name]))
            G.add_edges_from(itertools.product([reaction_name], product_list))        
        
        self.reaction_graph = G
        return G

    def isConsistent(self):
        """ Checks if the graph is consistent using the StoichiometryMatrix class. 

        Args:
            None. 

        Returns:
            True if consistent, False if inconsistent. 
        """
        sm_class = sm.StoichiometryMatrix(self.model)
        sm_class.buildMatrix()
        self.consistent = sm_class.isConsistent()

        return self.consistent

    def getInconsistentReactions(self):
        """ Find a set of inconsistent reactions.

        This method detects reactions from the cycles in a reaction graph
        under two criteria:

        1. Fina a cycle. If the cycle has a reaction with different number of
        reactants and products, choose reactions with the same number of
        reactants/products within the cycle. For example, if a cycle is
        A -> B + C // B -> D -> A, we choose 'D->A' as a problematic reaction. 

        2. Choose reactions with a self-loop (for example, A -> A + B)

        Args:
            None. 

        Returns:
            A set of possibly inconsistent reactions. 
        """
        reaction_list = [reaction.getId() for reaction in self.model.getListOfReactions()]
        if self.consistent == False:
            inconsistent_reactions = set()
            selfloop_reactions = set()
            for cycle in list(nx.simple_cycles(self.reaction_graph)):
                cycle_reaction = set(cycle).intersection(set(reaction_list))
                
                subG_list = [x for x in self.reaction_graph.edges if cycle_reaction.intersection(x) != set()]
                subG = nx.DiGraph()
                subG.add_edges_from(subG_list)

                sub_reaction_nodes = cycle_reaction
                sub_species_nodes = set(subG.nodes).difference(set(reaction_list))
                
                find_reactions = {reaction for reaction in sub_reaction_nodes \
                                         if subG.in_degree(reaction) - subG.out_degree(reaction) != 0}
                inconsistent_reactions = inconsistent_reactions.union(sub_reaction_nodes.difference(find_reactions))
                
                # find selfloop
                for tup in subG.edges:
                    if tuple(reversed(tup)) in subG.edges:
                        selfloop_reactions = selfloop_reactions.union( set(tup).intersection(sub_reaction_nodes) )
                
            self.inconsistent_set = inconsistent_reactions.union(selfloop_reactions)

        return self.inconsistent_set





