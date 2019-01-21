# -*- coding: utf-8 -*-
"""Doctrings on mass_inequality_graph.py.

-> High level description: What is the problem? Motivation of this problem 
(relationship between species in chemical reactions)

The MassInequalityGraph class represents a directed graph class with multiple edges between species. 
Each reactant is connected to a product, where each edge has a 'inequality' attribute, 
which is one of '=', '<', or '>'. This represents relative inequality of masses between
reactants and products, and helps us identify reactions that are not compatible. 

Example:
    You should already have an SBML model loaded, which we call 'model'.

    >> import mass_inequality_graph as mig
    >> mig_class = mig.MassInequalityGraph(model)
    >> mig_class.buildMassInequalityGraph()
"""


import constants as cn
import print_model as pm
import reaction_graph as rg
import stoichiometry_matrix as sm

import itertools  
import networkx as nx
import os
import tellurium as te
import tesbml


class MassInequalityGraph():
    """Creates a mass inequaltiy graph. 
    
    -> Here, more of software organization:
    Data structure or high level information:
    -> only use '<', and don't need '>'



    This class takes one input argument, an SBML model, and creates a 
    reaction graph class. From there, it creates another directed graph,
    this time with multiple edges between nodes. There are a few methods that
    identifies mass-inconsistent reactions. 

    -->> don't use this 
    Attributes:
        model (SBML model): An SBML model. The argument needed to create a reaction graph object. 
        reaction_graph (networkx.DiGraph): A reaction graph (directed graph) between reactants and products.
        imbalance_set(set): A set of reactions that are mass-inconsistent. 
    -->> don't use this

    """    
    def __init__(self, model):
        self.model = model

        # -->> don't use name 'class'; for example, instead of rg_class, reaction_graph (despite self.reaction_graph means something else)
        rg_class = rg.ReactionGraph(model)
        rg_class.buildReactionGraph()
        self.reaction_graph = rg_class.reaction_graph
        self.mass_inequality_graph = None # -> add explanation of it
        self.imbalance_set = None # -> add explanation of it

        # if not isinstance(self.model, tesbml.libsedml.Model):
        if self.model == None:
            raise TypeError("Model doesn't exist")

    def buildMassInequalityGraph(self):
        """Creates a mass inequality graph from a reaction graph. 
        
        Each reaction within the reaction graph has at least one reactant and one product.
        Here, we remove each reaction name node and directly connect reactants and products.
        At the same time, each edge has one of the two inequality attributes: '=' and '<'.
        Since the same edge between two nodes can have multiple attributes, we create a MultiDiGraph
        which will allow multiple edges between same pair of nodes. 
       
        Args:
            None. 

        Returns:
            mass_inequality_graph (networkx.classes.digraph.MultiDiGraph)    
        """
        reaction_list = [reaction.getId() for reaction in self.model.getListOfReactions()]
        MIG = nx.MultiDiGraph()
        for reaction in reaction_list:
            
            # skip the iteration if either reactant or product is empty
            # 'EmptySet' will still be included, such as in model #6, but it will be removed in MESGraph
            if (self.model.getReaction(reaction).getReactant(0) is None) | \
                (self.model.getReaction(reaction).getProduct(0) is None):
                continue
            
            in_nodes = [edge[0] for edge in self.reaction_graph.in_edges(reaction)]
            out_nodes = [edge[1] for edge in self.reaction_graph.out_edges(reaction)]
            
            # first, need to check if number of reactants/products is 1
            # next, need to confirm the stoichiometry is exactly 1.0
            if (len(in_nodes)==1) & (len(out_nodes)==1):
                if (self.model.getReaction(reaction).getReactant(0).getStoichiometry() == 1.0) & \
                    (self.model.getReaction(reaction).getProduct(0).getStoichiometry() == 1.0):
                    MIG.add_edge(in_nodes[0], out_nodes[0], inequality=cn.EQUAL, reaction=reaction)
                    MIG.add_edge(out_nodes[0], in_nodes[0], inequality=cn.EQUAL, reaction=reaction)      
            elif (len(in_nodes)==1) & (len(out_nodes)>1): 
                if self.model.getReaction(reaction).getReactant(0).getStoichiometry() == 1.0:
                    MIG.add_edges_from(itertools.product(in_nodes, out_nodes), inequality=cn.GREATERTHAN, \
                           reaction=reaction)
            elif (len(in_nodes)>1) & (len(out_nodes)==1): 
                if self.model.getReaction(reaction).getProduct(0).getStoichiometry() == 1.0:
                    MIG.add_edges_from(itertools.product(in_nodes, out_nodes), inequality=cn.LESSTHAN, \
                           reaction=reaction)

        self.mass_inequality_graph = MIG
        return MIG


    ### ->> come up with different name for this
    def buildMIG2(self):
        """Creates a MIG2, a maximally ordered graph by equivalence relations
        
        Each node of MIG2 is a set whoose elements have equal masses according to 
        the original mass inequality graph. The edges are created only if the species
        are strictlly less or great than each other. If unequal species show up within
        the same node, this reaction violates ordering and will be flagged. 
       
        Args:
            None. 

        Returns:
            MIG2 (networkx.classes.digraph.MultiDiGraph)    
        """

        ## -->>> use 
        reaction_list = [reaction.getId() for reaction in self.model.getListOfReactions()]
        ## -->>> use 'reactions' instead

        species_set = set(species.getId() for species in self.model.getListOfSpecies())
        edge_with_ineq = list(self.mass_inequality_graph.edges(data='inequality'))

        # choose reactions without boundary species
        reduced_edge_with_ineq = {edge for edge in edge_with_ineq if len(set(edge).intersection(species_set)) == 2}
        equal_edges = {(a,b,c) for a, b, c in reduced_edge_with_ineq if c == cn.EQUAL}
        unequal_edges = {(a,b,c) for a, b, c in reduced_edge_with_ineq if c != cn.EQUAL}
        model_species = species_set.intersection(set(self.mass_inequality_graph.nodes))

        equal_node_list = []

        num_equal_edges = len(equal_edges)
        remove_edges = equal_edges.copy()
        for edge in list(equal_edges):

            # take a (random) edge from the remove_edges set
            is_included = [(edge[0] in node_set) | (edge[1] in node_set) for node_set in equal_node_list]
            idx_included = [idx for idx, included in enumerate(is_included) if included]
            idx_not_included = [idx for idx, included in enumerate(is_included) if not included]
            
            # if there is no matching set, we create a new one
            if sum(is_included)==0:
                equal_node_list.append(frozenset({edge[0], edge[1]}))
            
            # if there is exactly one matching subset, we added this element to the set
            elif sum(is_included)==1:
                idx = idx_included[0]
                equal_node_list[idx] = equal_node_list[idx].union(frozenset({edge[0], edge[1]}))
                
            # if there are more than one matching subsets, we combine them and add the element
            else:
                merged_set = frozenset({edge[0], edge[1]})
                for idx in idx_included:
                    merged_set = merged_set.union(equal_node_list[idx])
                equal_node_list = [equal_node_list[idx] for idx in idx_not_included] + [merged_set]
            
            # remove an edge at every step
            remove_edges.remove(edge) 

        # identify model species that are not included in the equal_node_list. Create them as one-element sets?
        remaining_species_set = model_species.copy()
        for equal_node in equal_node_list:
            remaining_species_set = remaining_species_set.difference(equal_node)
        equal_node_list = equal_node_list + [frozenset([species]) for species in remaining_species_set]

        # now, create MIG2
        MIG2 = nx.MultiDiGraph()
        MIG2.add_nodes_from(equal_node_list)

        for edge in list(unequal_edges):
            is_icld_rct = [edge[0] in node_set for node_set in equal_node_list]
            is_icld_pdt = [edge[1] in node_set for node_set in equal_node_list]

            if is_icld_rct.index(True)==is_icld_pdt.index(True):
                print("Edge ", edge, " is inconsistent. ")
                print(edge[0], edge[2], edge[1], "and ", edge[0], cn.EQUAL, edge[1])
                print("")
            else:
                print("Edge", edge, " is applied to MIG2")
                MIG2.add_edge(equal_node_list[is_icld_rct.index(True)], 
                             equal_node_list[is_icld_pdt.index(True)])
        self.MIG2 = MIG2
        return MIG2


    # def findImbalancedReactions(self):
    #     """ Find a set of reactions that are inconsistent.

    #     This method identifies edges that have incompatible inequality attributes.
    #     For example, if two species 'A' and 'B' have two edges, each of which 
    #     has '=' and '<' respectively, we consider the relationship is mass-imbalanced. 

    #     Args:
    #         None. 

    #     Returns:
    #        imbalance_set (dictionary): A dictionary with species pairs and inconsistet attributes. 
    #     """
    #     imbalance_set = {}
    #     edge_list = list(set(self.mass_inequality_graph.edges(keys=False)))               
    #     multi_edge_list = [edge for edge in edge_list if self.mass_inequality_graph.number_of_edges(*edge)>1]

    #     for edge in multi_edge_list:
    #         edge_data = self.mass_inequality_graph.get_edge_data(*edge)
    #         ineq_set = {edge_data[key][cn.INEQUALITY] for key in edge_data}
    #         # length>1 means there are inconsistent relatioship between two species
    #         if len(ineq_set)>1:
    #             imbalance_pairs = []
    #             for ineq in list(ineq_set):
    #                 imbalance_pairs.append(edge[0] + ineq + edge[1])
    #             imbalance_set[edge[0]+'<->'+edge[1]] = imbalance_pairs

    #     self.imbalance_set = imbalance_set
    #     return imbalance_set





