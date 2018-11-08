# -*- coding: utf-8 -*-
"""Doctrings on graph_checker.py.

This script defines a graph checker class, which creates a set of undirected graphs
between reactants, products, and members of kinetic laws and assignment rules in
each reaction in a model. There are two types of graphs: MassFlowGraph and FullGraph.

MassFlowGraph connects only reactants and products in each reaction,
while FullGraph adds additional edges between products and species from
kinetic laws/assignment rules to the MassFlowGraph. 

Example:
    So far it is best to run on python mode or Jupyter Notebook/Lab if you want to test it. 
    You should already have an SBML model, which we call 'model' here.

        >> import graph_checker as gc
        >> onegraph = gc.GraphChecker(model)
        >> onegraph.isConnected() 
        >> onegraph.isMassConnected()

Attributes: /// now done yet
    module_level_variable1 (int): Module level variables may be documented in
        either the ``Attributes`` section of the module docstring, or in an
        inline docstring immediately following the variable.

        Either form is acceptable, but the two should not be mixed. Choose
        one convention to document module level variables and be consistent
        with it.

Todo:
    * Need to discuss if species from kinetic laws and assignments 
      (would you call them regulators?) 
      should be linked to products or reactants, or both.
    * Attributes?
    * Need to explain constants here? 
    * Discuss unittest for this class.
    * We may not need MassFlowGraph so need to discuss. 

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
STATUS = ["MODEL DOESN'T EXIST", "NO EDGES", "STRANDED GRAPH EXISTS", "CONNECTED"]
"""STATUS is a list of connectivity status of a model.

MODEL DOESN'T EXIST: Input argument (model) is None.
NO EDGES: Species are not connected to each other. Need to examine further.
STRANDED GRAPH EXISTS: One or more graphs are disconnected to another. Need to examine further.
CONNECTED: Model is normal in terms of connectivity. No need to check further. 
"""

# current_dir = os.getcwd()
# data_dir = os.path.abspath(os.path.join(current_dir, os.pardir, 'curated_data'))
# Trying to use the format in sphinx:http://www.sphinx-doc.org/en/master/

class GraphChecker():
    """Creates an undirected graph and checks connectivity.

    This class takes one input argument, an SBML model, and creates
    either a MassFlowGraph or FullGraph. The user will can obtain a graph
    by using getMassFlowGraph() or getFullGraph(), or just check 
    its connectivity by using isMassConnected() or isConnected(). 

    Attributes:
        model (SBML model): An SBML model. The argument needed to create a GraphChecker object. 
        graph (Graph object): A NetworkX Graph object.
        connected (str): Connectivity status of a model.

    """    
    def __init__(self, model):
        self.model = model

    def getMassFlowGraph(self):
        """Creates a MassFlowGraph and returns it.

        Note:
            If model is None, returns None. 

        Args:
            No specific arguments. 

        Returns:
            None if model doesn't exist, and a Graph object if it exists. 

        """
        if self.model == None:
            return None

        else:
            graph_set = set()
            for reaction_idx in range(self.model.getNumReactions()):
                reaction = self.model.getReaction(reaction_idx)
        
            reactant_set = set()
            for reactant_idx in range(reaction.getNumReactants()):
                species = reaction.getReactant(reactant_idx).getSpecies()
                if species != EMPTYSET:
                    reactant_set.add(species)

            product_set = set()
            for product_idx in range(reaction.getNumProducts()):
                species = reaction.getProduct(product_idx).getSpecies()
                if species != EMPTYSET:
                    product_set.add(species)

            graph_set = graph_set.union(itertools.product(reactant_set, product_set))

            G = nx.Graph()
            G.add_edges_from(graph_set)

            return G

    def isMassConnected(self):
        """Checks connectivity of MassFlowGraph and returns result.

        Args:
            No specific arguments. 

        Returns:
            One of members(str) from the STATUS list. 

        """
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
            
    def getFullGraph(self):
        """Creates a FullGraph and returns it.

        Note:
            If model is None, returns None. 

        Args:
            No specific arguments. 

        Returns:
            None if model doesn't exist, and a Graph object if it exists. 

        """
        if self.model == None:
            return None

        else:
            graph_edges = set()
            species_set = set()
            rules_set = set()

            # create a set of species
            for spec in self.model.getListOfSpecies():
                species_set.add(spec.getId())

            # create a set of assignment rules IDs
            for rule in self.model.getListOfRules():
                rules_set.add(rule.getId())
    
    
            for reaction_idx in range(self.model.getNumReactions()):
                reaction = self.model.getReaction(reaction_idx)
    
                reactants_set = set()
                for reactant_idx in range(reaction.getNumReactants()):
                    species = reaction.getReactant(reactant_idx).getSpecies()
                    if species != EMPTYSET:
                        reactants_set.add(species)

                products_set = set()
                for product_idx in range(reaction.getNumProducts()):
                    species = reaction.getProduct(product_idx).getSpecies()
                    if species != EMPTYSET:
                        products_set.add(species)

                # add regulators from Kinetic Law to connect with products
                kinetics = self.model.getReaction(reaction_idx).getKineticLaw()
                kinetic_formula = kinetics.getFormula()
                kinetic_atomic = set(REGEX.findall(kinetic_formula))   
                kinetic_atomic_subset = kinetic_atomic.intersection(species_set)
    
                graph_edges = graph_edges.union(itertools.product(kinetic_atomic_subset, products_set))
    
                # add regulators from Assignment Rules
                rules_subset = kinetic_atomic.intersection(rules_set)
                if bool(rules_subset):
                    for sub_rule in rules_subset:
                        rule_formula = self.model.getRule(sub_rule).getFormula()
                        rule_atomic = set(REGEX.findall(rule_formula))
                        rule_atomic_subset = rule_atomic.intersection(species_set)
                        graph_edges = graph_edges.union(itertools.product(rule_atomic_subset, products_set))
    

                graph_edges = graph_edges.union(itertools.product(reactants_set, products_set))

            G = nx.Graph()
            G.add_edges_from(graph_edges)

            return G

    def isConnected(self):
        """Checks connectivity of FullGraph and returns result.

        Args:
            No specific arguments. 

        Returns:
            One of members(str) from the STATUS list. 

        """
        self.graph = self.getFullGraph()

        if self.graph == None:
            self.connected = STATUS[0]
        elif nx.number_of_nodes(self.graph) == 0:
            self.connected = STATUS[1]
        elif nx.is_connected(self.graph) == False:
            self.connected = STATUS[2]
        else:
            self.connected = STATUS[3]
        
        return self.connected       


def display_statistics():
    
    """NEEDS TO BE RE-WRITTEN TO ADJUST TO THE GRAPHCHECKER CLASS.

    Runs graph_checker from 1 to 706 (except 596)
    and displays results and plots a bar graph
    
    Parameters:
    -----------
    None
    
    Returns:
    -----------
    null_graphs:disconnected_graphs, and wrong_format
    
    Displays a plot
    """
    
    wrong_format = []
    null_graphs = []
    disconnected_graphs = []
    run_list = set(range(1, 707)).difference({596})

    for i in run_list:
        graph_checker(i, wrong_format, null_graphs, disconnected_graphs)

    print("Number of null graphs: ", len(null_graphs))
    print("Number of disconnected graphs: ", len(disconnected_graphs))
    print("Finally, number of wrong formats(no models retrieved): ", len(wrong_format))
    print(wrong_format)
    
    bar_val = [707-len(null_graphs)-len(disconnected_graphs)-len(wrong_format),
              len(null_graphs), len(disconnected_graphs), len(wrong_format)]
    x_cat = ["Correct Models(" + str(bar_val[0]) + ")", 
             "No Graph(" + str(bar_val[1]) + ")",
             "Disconnected Graphs(" + str(bar_val[2]) + ")",
             "Couldn't Load(" + str(bar_val[3]) + ")"]
    plt.figure(figsize=(10,5))
    plt.title("Summary Statistics: Graph Checker")
    plt.bar(x_cat, bar_val, width=0.5)
    plt.show()
