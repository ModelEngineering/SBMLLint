# -*- coding: utf-8 -*-
"""Doctrings on bipartite_flow_graph.py.

The BipartiteFlowGraph class represents a directed graph between reactants and products
in each reaction. Each reactant is linked to a reaction object (a list of reaction name)
which is then connected to each product. 

Example:
    You should already have an SBML model loaded, which we call 'model'.

    >> import bipartite_flow_graph as fg
    >> sbml_graph = fg.BipartiteFlowGraph(model)

"""


import constants as cn

import itertools  
import networkx as nx
import os
import tellurium as te
import tesbml


class BipartiteFlowGraph():
    """Creates a directed mass flow graph.

    This class takes one input argument, an SBML model, and creates
    a mass flow graph called bipartite_flow_graph. The user can use a method
    to check its connectedness. 

    Attributes:
        model (SBML model): An SBML model. The argument needed to create a ModelGraph object. 
        connected (str): Connectivity status of a model.

        bipartite_flow_graph: A directed graph between reactants-reactions name-products 

    """    
    def __init__(self, model):
        self.model = model
        self.connected = None
        self.bipartite_flow_graph = None

        # if not isinstance(self.model, tesbml.libsedml.Model):
        if self.model == None:
            raise TypeError("Model doesn't exist")

    def buildBipartiteFlowGraph(self):
        """Creates a full stoichiometry matrix from a model.
        
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
        
        self.bipartite_flow_graph = G
        return G

    def isConnected(self):
        """ Checks if the graph if the graph is connected. 

        Args:
            None. 

        Returns:
            True if connected, False if not connected
        
        self.graph = self.getMassFlowGraph()

        if self.graph == None:
            self.connected = STATUS[0]
        elif nx.number_of_nodes(self.graph) == 0:
            self.connected = STATUS[1]
        elif nx.is_connected(self.graph) == False:
            self.connected = STATUS[2]
        else:
            self.connected = STATUS[3]
        
        return self.connected
"""