"""Mass Equality Set Graph (MESGraph)."""

from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule, MoleculeStoichiometry
from SBMLLint.common.reaction import Reaction
from SBMLLint.games import SOM
from SBMLLint.common.simple_sbml import SimpleSBML

import itertools
import networkx as nx


class MESGraph(nx.DiGraph):
	"""
    The MESGraph class represents a collection of SOMs as nodes
    and their inequality relationships as edges (arcs). 
    Mass inequality between SOMs from reactions can help us
    detect their relationship.
    Type I Error occurs when we find inequality between two molecules
    in the same SOM, because each element in a SOM has the same weight.
    Type II Error implies there is cyclism between molecules, such as
    A < B < C < ... < A, which is physically impossible. 
    """
	# all = []

    def __init__(self, soms):
        """
        :param list-SOM soms:
        """
        super(MESGraph, self).__init__()
        self.add_nodes_from(soms)
        self.identifier = self.makeId()

    def __repr__(self):
        return self.identifier
    
    def makeId(self):
        identifier = ":"
        for som in list(self.nodes):
            identifier = identifier + som.identifier + "+"
        return identifier
    
    def processUniUniReaction(self, reaction):
        """
        :param Reaction reactions:
        """
        if reaction.category != cn.REACTION_1_1:
            pass
        else:
            reactant_som = SOM.findSOM(reaction.reactants[0].molecule)
            product_som = SOM.findSOM(reaction.products[0].molecule)
            if reactant_som == product_som:
                pass
            else:
                new_som = SOM.merge(reaction)
                self.remove_node(reactant_som)
                self.remove_node(product_som)
                self.add_node(new_som)
                return new_som
    
    def processUniMultiReaction(self, reaction):
        """
        :param Reaction reaction:
        """
        if (reaction.category != cn.REACTION_1_n):
            pass
        else:
            destination = [SOM.findSOM(reaction.reactants[0].molecule)]
            source = [SOM.findSOM(product.molecule) for product in reaction.products]
            self.addArc(source, destination, reaction)

    def processMultiUniReaction(self, reaction):
        """
        :param Reaction reaction:
        """
        if (reaction.category != cn.REACTION_n_1):
            pass
        else:
            destination = [reaction.products[0].molecule]
            source = [reactant.molecule for reactant in reaction.reactants]
            self.addArc(source, destination, reaction)
    
    def addArc(self, source, destination, reaction):
        arcs = itertools.product(source, destination)
        for arc in arcs:
            if not self.checkTypeOneError(arc, reaction):
                self.add_edge(SOM.findSOM(arc[0]), SOM.findSOM(arc[1]), reaction=reaction)
            else:
                continue
    
    def checkTypeOneError(self, arc, reaction=None):
        som1 = SOM.findSOM(arc[0])
        som2 = SOM.findSOM(arc[1])
        if som1 == som2:
            print("We have Type I Error...")
            print(arc[0], " and ", arc[1], " have the same weight.")
            print("However, reaction -", reaction.label, "- implies ", arc[0], " < ", arc[1])
            print()
            return True
        else:
            return False
    
    def analyze(self, reactions):
        """
        :param list-Reaction reactions:
        """
        uniuni = []
        unimulti = []
        multiuni = []
        multimulti = []
        for reaction in reactions:
            if reaction.category == cn.REACTION_1_1:
                uniuni.append(reaction)
            elif reaction.category == cn.REACTION_1_n:
                unimulti.append(reaction)
            elif reaction.category == cn.REACTION_n_1:
                multiuni.append(reaction)
            elif reaction.category == cn.REACTION_n_n:
                multimulti.append(reaction)
            
        for reaction in uniuni:
            self.processUniUniReaction(reaction)
        for reaction in unimulti:
            self.processUniMultiReaction(reaction)
        for reaction in multiuni:
            self.processMultiUniReaction(reaction)
        return 



