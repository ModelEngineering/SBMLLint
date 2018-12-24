# -*- coding: utf-8 -*-
"""Doctrings on mass_inequality_graph.py.

The MassInequalityGraph class represents a directed graph class with multiple edges between species. 
Each reactant is connected to a product, where each edge has a 'inequality' attribute, 
which is one of '=', '>', and '<'. This represents relative inequality of masses between
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

    This class takes one input argument, an SBML model, and creates a 
    reaction graph class. From there, it creates another directed graph,
    this time with multiple edges between nodes. There are a few methods that
    identifies mass-inconsistent reactions. 

    Attributes:
        model (SBML model): An SBML model. The argument needed to create a reaction graph object. 
        reaction_graph (networkx.DiGraph): A reaction graph (directed graph) between reactants and products.
        imbalance_set(set): A set of reactions that are mass-inconsistent. 

    """    
    def __init__(self, model):
        self.model = model
        rg_class = rg.ReactionGraph(model)
        rg_class.buildReactionGraph()
        self.reaction_graph = rg_class.reaction_graph
        self.mass_inequality_graph = None
        self.imbalance_set = None

        # if not isinstance(self.model, tesbml.libsedml.Model):
        if self.model == None:
            raise TypeError("Model doesn't exist")

    def buildMassInequalityGraph(self):
        """Creates a mass inequality graph from a reaction graph. 
        
        Each reaction within the reaction graph has at least one reactant and one product.
        Here, we remove each reaction name node and directly connect reactants and products.
        At the same time, each edge has one of the three inequality attributes: '=', '<', and '>'.
        Since the same edge between two nodes can have multiple attributes, we create a MultiDiGraph
        which will allow multiple edges between same pair of nodes. 
       
        Args:
            None. 

        Returns:
            networkx.classes.digraph.MultiDiGraph    
        """
        reaction_list = [reaction.getId() for reaction in self.model.getListOfReactions()]
        MIG = nx.MultiDiGraph()
        for reaction in reaction_list:
            in_nodes = [edge[0] for edge in self.reaction_graph.in_edges(reaction)]
            out_nodes = [edge[1] for edge in self.reaction_graph.out_edges(reaction)]
              
            if (len(in_nodes)==1) & (len(out_nodes)==1):
                MIG.add_edge(in_nodes[0], out_nodes[0], inequality='=', reaction=reaction)
                MIG.add_edge(out_nodes[0], in_nodes[0], inequality='=', reaction=reaction)      
            elif (len(in_nodes)==1) & (len(out_nodes)>1): 
                MIG.add_edges_from(itertools.product(in_nodes, out_nodes), inequality='>', \
                           reaction=reaction)
            elif (len(in_nodes)>1) & (len(out_nodes)==1): 
                MIG.add_edges_from(itertools.product(in_nodes, out_nodes), inequality='<', \
                           reaction=reaction)

        self.mass_inequality_graph = MIG
        return MIG



    def findImbalancedReactions(self):
        """ Find a set of reactions that are inconsistent.

        This method identifies edges that have incompatible inequality attributes.
        For example, if two species 'A' and 'B' have two edges, each of which 
        has '=' and '<' respectively, we consider the relationship is mass-imbalanced. 

        Args:
            None. 

        Returns:
           imbalance_set (dictionary): A dictionary with species pairs and inconsistet attributes. 
        """
        imbalance_set = {}
        edge_list = list(set(self.mass_inequality_graph.edges(keys=False)))               
        multi_edge_list = [edge for edge in edge_list if self.mass_inequality_graph.number_of_edges(*edge)>1]

        for edge in multi_edge_list:
            edge_data = self.mass_inequality_graph.get_edge_data(*edge)
            ineq_set = {edge_data[key][cn.INEQUALITY] for key in edge_data}
            # length>1 means there are inconsistent relatioship between two species
            if len(ineq_set)>1:
                imbalance_pairs = []
                for ineq in list(ineq_set):
                    imbalance_pairs.append(edge[0] + ineq + edge[1])
                imbalance_set[edge[0]+'<->'+edge[1]] = imbalance_pairs

        self.imbalance_set = imbalance_set
        return imbalance_set





