"""Mass Equality Set Graph (MESGraph)."""

from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule, MoleculeStoichiometry
from SBMLLint.common.reaction import Reaction
from SBMLLint.games.som import SOM
from SBMLLint.common.simple_sbml import SimpleSBML

import collections
import itertools
import networkx as nx


PathComponents = collections.namedtuple('PathComponents',
                                        'node1 node2 reactions')


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

  def __init__(self, simple=None):
    """
    :param SimpleSBML simple:
    """
    super(MESGraph, self).__init__()
    self.simple = simple
    self.soms = self.initializeSOMs(simple)
    self.add_nodes_from(self.soms)
    self.identifier = self.makeId()
    self.type_one_error = False
    self.type_two_error = False
    self.type_one_errors = []
    self.type_two_errors = []

  def __repr__(self):
    return self.identifier

  def initializeSOMs(self, simple):
    """
    Create a list of one-molecule SOMs
    :param SimpleSBML simple:
    :return list-SOM:
    """
    soms = []
    if type(simple) == SimpleSBML:
      for molecule in simple.molecules:
        if molecule.name == cn.EMPTYSET:
          continue
        else:
          soms.append(SOM({molecule}))
    return soms

  def makeId(self):
    """
    Construct an identifier for the graph.
    :return str:
    """
    identifier = ""
    if self.edges:
      for edge in self.edges:
        identifier = identifier + str(edge[0]) + cn.ARC_ARROW + str(edge[1]) + "\n"
    for key, node in enumerate(nx.isolates(self)):
      identifier = identifier + str(node)
      if key < (len(list(nx.isolates(self)))-1):
          identifier = identifier + cn.KINETICS_SEPARATOR
    # Return the identifier
    return identifier

  def getNode(self, molecule):
    """
    Find a node(SOM) containing the given molecule.
    If no such SOM exists, return False
    :param Molecule molecule:
    :return SOM/False:
    """
    for som in list(self.nodes):
      for mole in som.molecules:
        if mole.name == molecule.name:
          return som
    return False

  def processUniUniReaction(self, reaction):
    """
    Process a 1-1 reaction to merge nodes.
    If no need to merge, return None.
    :param Reaction reactions:
    """
    if reaction.category != cn.REACTION_1_1:
      pass
    else:
      reactant_som = self.getNode(reaction.reactants[0].molecule)
      product_som = self.getNode(reaction.products[0].molecule)
      if reactant_som == product_som:
        return None
      else:
        new_som = reactant_som.merge(product_som)
        new_som.reactions.add(reaction)
        # TODO: if there are edges, need to also check them
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
    if reaction.category != cn.REACTION_1_n:
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
    if reaction.category != cn.REACTION_n_1:
      pass
    else:
      destination = [reaction.products[0].molecule]
      source = [reactant.molecule for reactant in reaction.reactants]
      self.addArc(source, destination, reaction)
      self.identifier = self.makeId()

  def addArc(self, source, destination, reaction):
    """
    Add arcs (edges) using two molecule lists (source/destination).
    :param list-Molecule source:
    :param list-Molecule destination:
    """
    arcs = itertools.product(source, destination)
    for arc in arcs:
      if not self.checkTypeOneError(arc, reaction):
        arc_source = self.getNode(arc[0])
        arc_destination = self.getNode(arc[1])
        # if there is already a preious reaction,
        if self.has_edge(arc_source, arc_destination):
          reaction_label = self.get_edge_data(arc_source, arc_destination)[cn.REACTION]
          # if reaction.label is not already included in the attribute,
          if reaction.label not in set(reaction_label):
            reaction_label = reaction_label + [reaction.label]
        else:
          reaction_label = [reaction.label]
        self.add_edge(arc_source, arc_destination, reaction=reaction_label)
      else:
        continue

  def getSOMPath(self, som, mole1, mole2):
    """
    Create an undirected graph between
    two molecules within a SOM
    and find the shortest path
    :param SOM som:
    :param Molecule mole1:
    :param Molecule mole2:
    :return PathComponents som_path:
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
      som_path.append(PathComponents(node1=path[idx], 
                                     node2=path[idx+1],
                                     reactions=edge_reactions))
    return som_path

  def addTypeOneError(self, mole1, mole2, reaction):
    """
    Add Type I Error components.
    All components of resulting PathComponents are str
    :param Molecule mole1:
    :param Molecule mole2:
    :param Reaction reaction:
    :return bool flag:
    """
    flag = False
    for component in self.type_one_errors:
      if (component.node1==mole1.name) and (component.node2==mole2.name):
        new_component = PathComponents(node1=mole1.name, 
                                       node2=mole2.name,
                                       reactions=component.reactions+[reaction.label])
        self.type_one_errors.remove(component)
        self.type_one_errors.append(new_component)
        flag = True
        break
    if not flag:
      self.type_one_errors.append(PathComponents(node1=mole1.name, 
                                                 node2=mole2.name,
                                                 reactions=[reaction.label]))
      flag = True
    return flag

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
    som1 = self.getNode(arc[0])
    som2 = self.getNode(arc[1])
    if som1 == som2:
      # if error_details:
      #   print("We have a Type I Error...")
      #   print(arc[0], " and ", arc[1], " have the same weight.")
      #   som_path = self.getSOMPath(som1, arc[0], arc[1])
      #   for pat in som_path:
      #     print(pat.node1 + " = " + pat.node2 + " from ", end="")
      #     for r in pat.reactions:
      #       reaction = self.simple.getReaction(r)
      #       print(reaction.makeIdentifier(is_include_kinetics=False))
      #   print("\nHowever, the following reaction") 
      #   print(inequality_reaction.makeIdentifier(is_include_kinetics=False)) 
      #   print("implies ", arc[0], cn.LESSTHAN, arc[1])
      #   print()
      # add type I error to self.type_one_errors
      self.addTypeOneError(arc[0], arc[1], inequality_reaction)
      return True
    else:
      return False

  def checkTypeTwoError(self, error_details=True):
    """
    Check Type II Error (cycles) of a MESGraph.
    If there is at least one cycle, 
    report an error message, related reactions
    and return True.
    If there is no cycle, return False. 
    :return bool:
    """
    graph = nx.DiGraph()
    graph.add_edges_from(self.edges)
    cycles = list(nx.simple_cycles(graph))
    if len(cycles) == 0:
      return False
    else:
      if error_details:
        print("We have a Type II Error...\n")
        for cycle in cycles:
          cycle_nodes = []
          for node in cycle:
            cycle_nodes.append(node)
          for node in cycle_nodes:
            print(node, cn.LESSTHAN, end=" ")
          print(cycle_nodes[0], "\n")
          #
          for node_idx in range(0, len(cycle_nodes)-1):
            arc_source = cycle_nodes[node_idx]
            arc_destination = cycle_nodes[node_idx+1]
            print(arc_source, cn.LESSTHAN, arc_destination, " by")
            reaction_label = self.get_edge_data(arc_source, arc_destination)[cn.REACTION]
            for r_label in reaction_label:
              print(r_label + "\n")
          # last arc which completes the cycle
          arc_source = cycle_nodes[len(cycle_nodes)-1]
          arc_destination = cycle_nodes[0]
          print(arc_source, cn.LESSTHAN, arc_destination, " by")
          reaction_label = self.get_edge_data(arc_source, arc_destination)[cn.REACTION]
          for r_label in reaction_label:
            print(r_label + "\n")
        #
        if not self.type_two_error:
          self.type_two_error = True
        return True

  def analyze(self, reactions, error_details=True):
    """
    Sort list of reactions and process them.
    Add arcs or sending error messages using
    checkTypeOneError or checkTypeTwoError.
    :param list-Reaction reactions:
    """
    # Associate the reaction category with the function
    # that processes that category
    reaction_dic = {
        cn.REACTION_1_1: self.processUniUniReaction,
        cn.REACTION_1_n: self.processUniMultiReaction,
        cn.REACTION_n_1: self.processMultiUniReaction,
        }
    # Process each type of reaction
    for category in reaction_dic.keys():
      for reaction in [r for r in reactions if r.category == category]:
        func = reaction_dic[category]
        func(reaction)
    if error_details:
      # print("Type I Errors:", self.type_one_errors)
      for error_path in self.type_one_errors:
        #print(error_path.node1 + " and " + error_path.node2 + " have same weight.")
        node1 = self.simple.getMolecule(error_path.node1)
        node2 = self.simple.getMolecule(error_path.node2)
        som = self.getNode(node1)
        som_path = self.getSOMPath(som, node1, node2)
        for pat in som_path:
          print(pat.node1 + " = " + pat.node2 + " by reaction(s):")
          for r in pat.reactions:
            som_reaction = self.simple.getReaction(r)
            print(som_reaction.makeIdentifier(is_include_kinetics=False))
        print("\nHowever, the following reaction(s)") 
        for arc_reaction in error_path.reactions:
          print(self.simple.getReaction(arc_reaction).makeIdentifier(is_include_kinetics=False)) 
        print("imply " + node1.name, cn.LESSTHAN, node2.name)
        print("------------------------------------")

    self.checkTypeTwoError(error_details)
    #
    return self
