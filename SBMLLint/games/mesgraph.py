"""Mass Equality Set Graph (MESGraph)."""

from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule, MoleculeStoichiometry
from SBMLLint.common.reaction import Reaction
from SBMLLint.games.som import SOM
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
        identifier = ""
        if self.edges:
            for edge in self.edges:
                identifier = identifier + str(edge[0]) + cn.ARC_ARROW + str(edge[1]) + "\n"

        for key, node in enumerate(nx.isolates(self)):
            identifier = identifier + str(node)
            if key < (len(list(nx.isolates(self)))-1):
                identifier = identifier + cn.KINETICS_SEPARATOR
                
        return identifier
    
    def processUniUniReaction(self, reaction):
        """
        Process a 1-1 reaction to update nodes.
        Uses SOM.merge(reaction) method. 
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
                self.identifier = self.makeId()
                return new_som
    
    def processUniMultiReaction(self, reaction):
        """
        Process a 1-n reaction to add arcs.
        Since the mass of reactant is greater than
        that of each product, it adds arcs by
        addArc(source=products, destination=reactant). 
        :param Reaction reaction:
        """
        if (reaction.category != cn.REACTION_1_n):
            pass
        else:
            destination = [reaction.reactants[0].molecule]
            source = [product.molecule for product in reaction.products]
            self.addArc(source, destination, reaction)
            self.identifier = self.makeId()

    def processMultiUniReaction(self, reaction):
        """
        Process a n-1 reaction to add arcs.
        Since the mass of product is greater than
        that of each reactant, it adds arcs by
        addArc(source=reactants, destination=product). 
        :param Reaction reaction:
        """
        if (reaction.category != cn.REACTION_n_1):
            pass
        else:
            destination = [reaction.products[0].molecule]
            source = [reactant.molecule for reactant in reaction.reactants]
            self.addArc(source, destination, reaction)
            self.identifier = self.makeId()
    
    def addArc(self, source, destination, reaction):
        """
        Check Type I Error and Add arcs 
        using two list of molecules (source/destination).
        :param list-Molecule source:
        :param list-Molecule destination:
        """
        arcs = itertools.product(source, destination)
        for arc in arcs:
            if not self.checkTypeOneError(arc, reaction):
                self.add_edge(SOM.findSOM(arc[0]), SOM.findSOM(arc[1]), reaction=reaction)
            else:
                continue
    
    def checkTypeOneError(self, arc, inequality_reaction=None):
        """
        Check Type I Error of an arc.
        If both source and destination are found
        in the same SOM, send error message and return True.
        If not, return False.
        :param tuple-Molecule arc:
        :param Reaction inequality_reaction:
        :return bool:
        """
        som1 = SOM.findSOM(arc[0])
        som2 = SOM.findSOM(arc[1])
        if som1 == som2:
            print("We have Type I Error...")
            print(arc[0], " and ", arc[1], " have the same weight by")
            for equality_reaction in list(som1.reactions):
                print(equality_reaction)
            ##
            print("\nHowever, reaction \"", inequality_reaction, 
                  "\" implies ", arc[0], " < ", arc[1])
            print()
            return True
        else:
            return False
    
    def analyze(self, reactions):
        """
        Sort list of reactions and process them.
        Add arcs or sending error messages using
        checkTypeOneError or checkTypeTwoError.
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
            
        return self
