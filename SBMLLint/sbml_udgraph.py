# -*- coding: utf-8 -*-
"""Doctrings on sbml_udgraph.py

This script defines a model graph class, which creates a set of undirected graphs
between reactants, products, and members of kinetic laws and assignment rules in
each reaction in a model. There are two types of graphs: mass_graph and full_graph.

mass_graph connects only reactants and products in each reaction,
while full_graph adds additional edges between products and species from
kinetic laws/assignment rules to the mass_graph

Example:
    So far it is best to run on python mode or Jupyter Notebook/Lab if you want to test it. 
    You should already have an SBML model, which we call 'model' here.

        >> import sbml_udgraph as su
        >> onegraph = us.ModelGraph(model)
        >> onegraph.isConnected() 
        >> onegraph.isMassConnected()

"""

import itertools  
import os

import matplotlib.pyplot as plt
import networkx as nx
import tellurium as te
import re
import roadrunner
import sympy
import sys
import tesbml


EMPTYSET = 'EmptySet'
PATTERN = r'[A-Z0-9a-z_]{1,}'
REGEX = re.compile(PATTERN)
STATUS = ["NULL GRAPH", "DISCONNECTED", "CONNECTED"]
"""STATUS is a list of connectivity status of a model.

Null Graph: No species is connected to each other. Need to examine further.
Disconnected: One or more graphs are stranded. Need to examine further.
CONNECTED: Model is normal in terms of connectivity. No need to check further. 
"""

class ModelGraph:
    """Creates an undirected graph between reactants, products, and regulators
    and checks connectivity.

    This class takes one input argument, an SBML model, and creates
    a mass flow graph and full graph. The user can check its connectivity by
    either isMassConnected() or isConnected(). 

    Attributes:
        model (SBML model): An SBML model. The argument needed to create a ModelGraph object. 
        connected (str): Connectivity status of a model.

        mass_graph: An undirected graph created only with reactant-product pairs.

        full_graph: An undirected graph created with reactant-product and regulator-product pairs.  

    """    
    def __init__(self, model):
        self.model = model

        if self.model == None:
            raise TypeError("Model doesn't exist")

        mass_graph_edges = set()
        full_graph_edges = set()

        # species and rule sets are used to construct full graph
        species_set = {spec.getId() for spec in self.model.getListOfSpecies()}
        rule_set = {rule.getId() for rule in self.model.getListOfRules()}

        for reaction_idx in range(self.model.getNumReactions()):
            reaction = self.model.getReaction(reaction_idx) 

            reactant_set = {species.getSpecies() for species in reaction.getListOfReactants() if species.getSpecies() != EMPTYSET}
            product_set = {species.getSpecies() for species in reaction.getListOfProducts() if species.getSpecies() != EMPTYSET}

            # collect edges for mass graph
            mass_graph_edges = mass_graph_edges.union(itertools.product(reactant_set, product_set))

            kinetics = self.model.getReaction(reaction_idx).getKineticLaw()
            kinetic_formula = kinetics.getFormula()
            kinetic_atomic = set(REGEX.findall(kinetic_formula))   
            kinetic_atomic_subset = kinetic_atomic.intersection(species_set).difference(reactant_set).difference(product_set)

            full_graph_edges = full_graph_edges.union(mass_graph_edges)
            full_graph_edges = full_graph_edges.union(itertools.product(kinetic_atomic_subset, product_set))
    
            # add regulators from Assignment Rules
            rule_subset = kinetic_atomic.intersection(rule_set)
            if bool(rule_subset):
                for sub_rule in rule_subset:
                    rule_formula = self.model.getRule(sub_rule).getFormula()
                    rule_atomic = set(REGEX.findall(rule_formula))
                    rule_atomic_subset = rule_atomic.intersection(species_set).difference(reactant_set).difference(product_set)

                    full_graph_edges = full_graph_edges.union(itertools.product(rule_atomic_subset, product_set))
    
            # collect edges for full graph
            full_graph_edges = full_graph_edges.union(itertools.product(reactant_set, product_set))

        mass_graph = nx.Graph()
        full_graph = nx.Graph()

        mass_graph.add_edges_from(mass_graph_edges)
        full_graph.add_edges_from(full_graph_edges)

        self.mass_graph = mass_graph
        self.full_graph = full_graph


    def isMassConnected(self):
        """Checks connectivity of a mass graph and returns result.

        Args:
            No specific arguments. 

        Returns:
            One of members(str) from the STATUS list. 

        """

        if nx.number_of_nodes(self.mass_graph) == 0:
            self.mass_connected = STATUS[0]
        elif nx.is_connected(self.mass_graph) == False:
            self.mass_connected = STATUS[1]
        else:
            self.mass_connected = STATUS[2]

        return self.mass_connected
            

    def isConnected(self):
        """Checks connectivity of a full graph and returns result.

        Args:
            No specific arguments. 

        Returns:
            One of members(str) from the STATUS list. 

        """

        if nx.number_of_nodes(self.full_graph) == 0:
            self.full_connected = STATUS[0]
        elif nx.is_connected(self.full_graph) == False:
            self.full_connected = STATUS[1]
        else:
            self.full_connected = STATUS[2]

        return self.full_connected
