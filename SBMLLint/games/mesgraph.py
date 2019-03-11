"""Mass Equality Set Graph (MESGraph)."""

from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule, MoleculeStoichiometry
from SBMLLint.common.reaction import Reaction
from SBMLLint.games.som import SOM
from SBMLLint.common.simple_sbml import SimpleSBML

import collections
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

  def __init__(self, simple=None):
    """
    :param SimpleSBML simple:
    """
    super(MESGraph, self).__init__()
    self.simple = simple
    self.soms = self.initializeSOMs(simple)
    self.add_nodes_from(self.soms)
    self.identifier = self.makeId()
    self.multimulti_reactions = []
    self.type_one_error = False
    self.type_two_error = False
    self.type_one_errors = []
    self.type_two_errors = []
    self.type_three_errors = []

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

  def mergeNodes(self, som1, som2, reaction):
    """
    Merge two nodes (SOMs).
    Update arcs if applicable. 
    :param SOM som1:
    :param SOM som2:
    :param Reaction reaction:
    :return SOM new_som:
    """
    new_som = som1.merge(som2)
    new_som.reactions.add(reaction)
    for som in [som1, som2]:
      for edge in list(self.in_edges(som)):
        remaining_som = edge[0]
        reaction_label = self.get_edge_data(edge[0], edge[1])[cn.REACTION]
        self.add_edge(remaining_som, new_som, reaction=reaction_label)  
      for edge in list(self.out_edges(som)):
        remaining_som = edge[1]
        reaction_label = self.get_edge_data(edge[0], edge[1])[cn.REACTION]
        self.add_edge(new_som, remaining_som, reaction=reaction_label) 
    self.remove_nodes_from([som1, som2])
    if not self.has_node(new_som):
      self.add_node(new_som)
    return new_som

  def processUniUniReaction(self, reaction):
    """
    Process a 1-1 reaction to merge nodes.
    If no need to merge, return None.
    :param Reaction reaction:
    """
    if reaction.category != cn.REACTION_1_1:
      pass
    else:
      reactant_som = self.getNode(reaction.reactants[0].molecule)
      product_som = self.getNode(reaction.products[0].molecule)
      if reactant_som == product_som:
        return None
      else:
        # new_som = reactant_som.merge(product_som)
        # new_som.reactions.add(reaction)
        # # TODO: if there are edges, need to also check them: self.mergeNodes(som1, som2, reaction)
        # self.remove_node(reactant_som)
        # self.remove_node(product_som)
        # self.add_node(new_som)
        new_som = self.mergeNodes(reactant_som, product_som, reaction)
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
      arcs = itertools.product(source, destination)
      for arc in arcs:
        if not self.checkTypeOneError(arc, reaction):
          som_source = self.getNode(arc[0])
          som_destination = self.getNode(arc[1])
          self.addArc(som_source, som_destination, reaction)
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
      arcs = itertools.product(source, destination)
      for arc in arcs:
        if not self.checkTypeOneError(arc, reaction):
          som_source = self.getNode(arc[0])
          som_destination = self.getNode(arc[1])
          self.addArc(som_source, som_destination, reaction)
      self.identifier = self.makeId()

  def addMultiMultiReaction(self, reaction=None):
    """
    Add a multi-multi reaction to self.multimulti_reactions
    :param reaction Reaction:
    :return bool:
    """
    if reaction not in self.multimulti_reactions:
      self.multimulti_reactions.append(reaction)

  def addTypeThreeError(self, som1, som2, reaction):
    """
    Add Type III Error components to self.type_three_errors
    All components of resulting PathComponents are str
    :param SOM som1:
    :param SOM som2:
    :param Reaction reaction:
    :return bool flag:
    """
    flag = False
    for component in self.type_three_errors:
      if (component.node1==som1) and (component.node2==som2):
        new_component = cn.PathComponents(node1=som1, 
                                          node2=som2,
                                          reactions=component.reactions+[reaction.label])
        self.type_three_errors.remove(component)
        self.type_three_errors.append(new_component)
        flag = True
        break
    if not flag:
      self.type_three_errors.append(cn.PathComponents(node1=som1, 
                                                      node2=som2,
                                                      reactions=[reaction.label]))
      flag = True
    return flag  

  def checkTypeThreeError(self, som1, som2, reaction):
    """
    Check type III error, which is when
    we cannot merge two nodes because there is an arc.
    Add the error to type three error. 
    :param SOM som1:
    :param SOM som2:
    :param reaction Reaction:
    :return bool:
    """
    if self.has_edge(som1, som2):
      self.addTypeThreeError(som1, som2, reaction)
      return True
    elif self.has_edge(som2, som1):
      self.addTypeThreeError(som2, som1, reaction)
      return True
    else:
      return False

  # def checkTypeFourError(self, som1, som2, reaction):
  #   """
  #   Check type IV error, which is when
  #   we cannot add an arc because there is an arc.
  #   Add the error to type three error.
  #   :param SOM som1:
  #   :param SOM som2:
  #   :param reaction Reaction:
  #   :return bool:
  #   """
  #   if som1 == som2:
  #     self.addTypeFourError(som1, som2, reaction)
  #     return True
  #   return False

  def reduceReaction(self, reaction):
    """
    Reduce the given reaction 
    :param Reaction reaction:
    :return False/Reaction reaction:
    """
    if reaction.category != cn.REACTION_n_n:
        return False
    # Reduces the reaction by examining for each SOM
    for som in list(self.nodes):
      reactants_in = collections.deque([mole_stoich for mole_stoich in  
                           reaction.reactants if 
                           self.getNode(mole_stoich.molecule)==som])
      reactants_out = [mole_stoich for mole_stoich in  
                      reaction.reactants if 
                      self.getNode(mole_stoich.molecule)!=som]
      products_in = collections.deque([mole_stoich for mole_stoich in  
                          reaction.products if 
                          self.getNode(mole_stoich.molecule)==som])
      products_out = [mole_stoich for mole_stoich in  
                     reaction.products if 
                      self.getNode(mole_stoich.molecule)!=som]
      #
      while reactants_in and products_in:
        reactant = reactants_in[0]
        product = products_in[0]
        if reactant.stoichiometry > product.stoichiometry:
          reactants_in[0] = MoleculeStoichiometry(reactant.molecule,
                                                  reactant.stoichiometry - product.stoichiometry)
          products_in.popleft()
        elif reactant.stoichiometry < product.stoichiometry:
          products_in[0] = MoleculeStoichiometry(product.molecule,
                                                 product.stoichiometry - reactant.stoichiometry)
          reactants_in.popleft()
        else:
          reactants_in.popleft()
          products_in.popleft()
      reactants = list(reactants_in) + reactants_out
      products = list(products_in) + products_out
    #  
    if (len(reaction.reactants) > len(reactants)) | \
    (len(reaction.products) > len(products)):
      reduced = True
      reaction.reactants = reactants
      reaction.products = products
    #
    reaction.identifier = reaction.makeIdentifier()
    reaction.category = reaction._getCategory() 
    return reaction

  def processMultiMultiReaction(self, reaction):
    """
    Process a multi-multi reaction.
    :param Reaction reaction:
    :return bool flag:
    """
    # need to redude the reaction first

    reduced_reaction = reaction
    reactant_soms = list({self.getNode(ms.molecule) for ms in reduced_reaction.reactants})
    product_soms = list({self.getNode(ms.molecule) for ms in reduced_reaction.products})
    if (len(reactant_soms)>1) and (len(product_soms)>1):
      return False
    #
    if (len(reactant_soms)==1) and (len(product_soms)==1):
      reactant_stoichiometry = sum([ms.stoichiometry for ms in reaction.reactants])
      product_stoichiometry = sum([ms.stoichiometry for ms in reaction.products])
      reactant_som = reactant_soms[0]
      product_som = product_soms[0]
      if reactant_stoichiometry == product_stoichiometry:
        if reactant_som != product_som:
          if not self.checkTypeThreeError(reactant_som, product_som, reaction):
            self.mergeNodes(reactant_som, product_som, reaction)
      # Add reactant_som -> product_som
      elif reactant_stoichiometry > product_stoichiometry:
        self.addArc(reactant_som, product_som, reaction)
      # Add product_som -> reactant_som
      else:
        self.addArc(product_som, reactant_som, reaction)
      return True
    else:
      if (len(reactant_soms)==1) and (len(product_soms)>1):
        som_arcs = itertools.product(product_soms, reactant_soms)
      elif (len(reactant_soms)>1) and (len(product_soms)==1):
        som_arcs = itertools.product(reactant_soms, product_soms)
      for arc in som_arcs:
        self.addArc(arc[0], arc[1], reaction)
      return True
    # return false if the reaction was not processed
    return False

  def addArc(self, arc_source, arc_destination, reaction):
    """
    Add a single arc (edge) using two SOMs and reaction.
    :param SOM arc_source:
    :param SOM arc_destination:
    :param Reaction reaction:
    """
    # if there is already a preious reaction,
    if self.has_edge(arc_source, arc_destination):
      reaction_label = self.get_edge_data(arc_source, arc_destination)[cn.REACTION]
      # if reaction.label is not already included in the attribute,
      if reaction.label not in set(reaction_label):
        reaction_label = reaction_label + [reaction.label]
    else:
      reaction_label = [reaction.label]
    # overwrite the edge with new reactions set
    self.add_edge(arc_source, arc_destination, reaction=reaction_label)

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
      som_path.append(cn.PathComponents(node1=path[idx], 
                                     node2=path[idx+1],
                                     reactions=edge_reactions))
    return som_path

  def printSOMPath(self, molecule_name1, molecule_name2):
    """
    Print out shortest SOM path between two molecules.
    Arguments are str and both molecules sholud be in the 
    same SOM.
    :param str molecule_name1:
    :param str molecule_name2:
    :return bool
    """
    som1 = self.getNode(self.simple.getMolecule(molecule_name1))
    som2 = self.getNode(self.simple.getMolecule(molecule_name2))
    if som1 != som2:
      return False
    else:
      # add case when molecule_name1 == molecule_name2
      if molecule_name1 == molecule_name2:
        print("Clearly,", molecule_name1, cn.EQUAL, molecule_name2)
      else:
        som_path = self.getSOMPath(som1, 
                                   self.simple.getMolecule(molecule_name1), 
                                   self.simple.getMolecule(molecule_name2))
        for pat in som_path:
          print("\n" + pat.node1, cn.EQUAL, pat.node2 + " by reaction(s):")
          for r in pat.reactions:
            som_reaction = self.simple.getReaction(r)
            print(som_reaction.makeIdentifier(is_include_kinetics=False))
      return True

  def addTypeOneError(self, mole1, mole2, reaction):
    """
    Add Type I Error components to self.type_one_errors
    All components of resulting PathComponents are str
    :param Molecule mole1:
    :param Molecule mole2:
    :param Reaction reaction:
    :return bool flag:
    """
    flag = False
    for component in self.type_one_errors:
      if (component.node1==mole1.name) and (component.node2==mole2.name):
        new_component = cn.PathComponents(node1=mole1.name, 
                                          node2=mole2.name,
                                          reactions=component.reactions+[reaction.label])
        self.type_one_errors.remove(component)
        self.type_one_errors.append(new_component)
        flag = True
        break
    if not flag:
      self.type_one_errors.append(cn.PathComponents(node1=mole1.name, 
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
      self.addTypeOneError(arc[0], arc[1], inequality_reaction)
      return True
    else:
      return False

  def addTypeTwoError(self, cycle):
    """
    Add Type II Error components to self.type_two_errors
    which is a list of lists
    All components of resulting PathComponents are str
    :param list-SOM cycle:
    """
    # exceptionally, here PathComponents are
    # node1=[], node2=[], reactions=[] and their index
    # of each component will match. All elements within nodes
    # are in the same SOM
    error_cycle = []
    for node_idx in range(len(cycle)-1):
      som1 = cycle[node_idx]
      som2 = cycle[node_idx+1]
      som1_moles = {mole.name for mole in list(som1.molecules)}
      som2_moles = {mole.name for mole in list(som2.molecules)}
      reactions = self.get_edge_data(som1, som2)[cn.REACTION]
      # all reactions (in an edge), should create a single PathComponent
      nodes1 = []
      nodes2 = []
      reaction_labels = []
      for r in reactions:
        reaction = self.simple.getReaction(r)
        if reaction.category == cn.REACTION_n_1:
          sources = {r.molecule.name for r in reaction.reactants}
          destinations = {p.molecule.name for p in reaction.products}
        elif reaction.category == cn.REACTION_1_n:
          sources = {p.molecule.name for p in reaction.products}
          destinations = {r.molecule.name for r in reaction.reactants}
        # for any reaction that addes arcs, len(nodes2)==1
        node2 = list(destinations.intersection(som2_moles))[0]
        for node1 in list(sources.intersection(som1_moles)):
          nodes1.append(node1)
          nodes2.append(node2)
          reaction_labels.append(reaction.label)
      error_cycle.append(cn.PathComponents(node1=nodes1, 
                                           node2=nodes2,
                                           reactions=reaction_labels))
    som1 = cycle[-1]
    som2 = cycle[0]
    som1_moles = {mole.name for mole in list(som1.molecules)}
    som2_moles = {mole.name for mole in list(som2.molecules)}
    reactions = self.get_edge_data(som1, som2)[cn.REACTION]
    # all reactions (in an edge), should create a single PathComponent
    nodes1 = []
    nodes2 = []
    reaction_labels = []
    for r in reactions:
      reaction = self.simple.getReaction(r)
      if reaction.category == cn.REACTION_n_1:
        sources = {r.molecule.name for r in reaction.reactants}
        destinations = {p.molecule.name for p in reaction.products}
      elif reaction.category == cn.REACTION_1_n:
        sources = {p.molecule.name for p in reaction.products}
        destinations = {r.molecule.name for r in reaction.reactants}
        # for any reaction that addes arcs, len(nodes2)==1
      node2 = list(destinations.intersection(som2_moles))[0]
      for node1 in list(sources.intersection(som1_moles)):
        nodes1.append(node1)
        nodes2.append(node2)
        reaction_labels.append(reaction.label)
    error_cycle.append(cn.PathComponents(node1=nodes1, 
                                        node2=nodes2,
                                        reactions=reaction_labels))
    self.type_two_errors.append(error_cycle)  

  def checkTypeTwoError(self):
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
      for cycle in cycles:
        self.addTypeTwoError(cycle)
        if not self.type_two_error:
          self.type_two_error = True
      return True

  def analyze(self, reactions=None, error_details=True):
    """
    Sort list of reactions and process them.
    Add arcs or sending error messages using
    checkTypeOneError or checkTypeTwoError.
    :param list-Reaction reactions:
    """
    if reactions is None:
      reactions = self.simple.reactions
    # Associate the reaction category with the function
    # that processes that category
    reaction_dic = {
        cn.REACTION_1_1: self.processUniUniReaction,
        cn.REACTION_1_n: self.processUniMultiReaction,
        cn.REACTION_n_1: self.processMultiUniReaction,
        cn.REACTION_n_n: self.addMultiMultiReaction,
        }
    # Process each type of reaction
    for category in reaction_dic.keys():
      for reaction in [r for r in reactions if r.category == category]:
        func = reaction_dic[category]
        func(reaction)
    #
    self.checkTypeTwoError()    
    #
    if error_details:
      if (len(self.type_one_errors)==0) and (len(self.type_two_errors)==0):
        print("No mass balance error found.")
      #
      for error_path in self.type_one_errors:
        self.printSOMPath(error_path.node1, error_path.node2)
        print("\nHowever, the following reaction(s)") 
        for arc_reaction in error_path.reactions:
          print(self.simple.getReaction(arc_reaction).makeIdentifier(is_include_kinetics=False)) 
        print("imply " + error_path.node1, cn.LESSTHAN,  error_path.node2)
        print("------------------------------------")
      print("************************************")
      #
      # print("We Do have type II Errors", self.type_two_errors)
      for cycle in self.type_two_errors:
        for idx, path_comp in enumerate(cycle):
          nodes1 = collections.deque(path_comp.node1)
          nodes2 = collections.deque(path_comp.node2)
          if idx < len(cycle)-1:
            next_nodes1 = collections.deque(cycle[idx+1].node1)
          else:
            next_nodes1 = collections.deque(cycle[0].node1)
          reactions = collections.deque(path_comp.reactions)
          # print SOM path between node elements
          if len(nodes1)>1:
            for node_idx in range(len(nodes1)-1):
              self.printSOMPath(nodes1[node_idx], nodes1[node_idx+1])
          if not set(nodes2).intersection(set(next_nodes1)):
            self.printSOMPath(nodes2[0], next_nodes1[0])
         #
          while nodes1:
            print()
            print(nodes1[0] + " < " + nodes2[0] + " by reaction:")
            arc_reaction = self.simple.getReaction(reactions[0])
            print(arc_reaction.makeIdentifier(is_include_kinetics=False))
            nodes1.popleft()
            nodes2.popleft()
            reactions.popleft()
        print("------------------------------------")
      print("************************************")
      #
    for multimulti in self.multimulti_reactions:
      self.processMultiMultiReaction(multimulti)
    if self.type_three_errors:
      print("We have type III errors\n", self.type_three_errors)
    else:
      print("Unfortunately(?), We don't have type III errors")
    #
    self.identifier = self.makeId()
    return self
