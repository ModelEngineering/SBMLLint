"""
Reporting class for GAMES Plus (GAMES_PP) algorithm
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule, MoleculeStoichiometry
from SBMLLint.common.reaction import Reaction
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.games.games_pp import SOMStoichiometry, SOMReaction, GAMES_PP, TOLERANCE
from SBMLLint.games.som import SOM
from SBMLLint.common import simple_sbml

import collections
import networkx as nx
import numpy as np
import pandas as pd

NULL_STR = ""
ReactionOperation = collections.namedtuple("ReactionOperation", 
    "reaction operation")
NUM_STAR = 70
PARAGRAPH_DIVIDER = "\n" + "-"*NUM_STAR + "\n"
REPORT_DIVIDER = "\n" + "*"*NUM_STAR + "\n"


class SimplifiedReaction(object):

  def __init__(self, reactants, products, label, mesgraph):
    """
    :param list-MoleculeStoichiometry reactants:
    :param list-MoleculeStoichiometry products:
    :param str label:
    """
    self.reactants = reactants
    self.products = products
    self.label = label
    self.identifier = self.makeIdentifier()
    self.mesgraph = mesgraph

  def __repr__(self):
    return self.identifier
    
  def makeIdentifier(self):
    """
    Provides a string representation of the reaction
    :param bool is_include_kinetics: include the kinetics formula
    :return str:
    """
    def makeStoichiometryString(molecule_stoichiometry):
      num = molecule_stoichiometry.stoichiometry
      if np.isclose(num, 1.0):
        return ''
      else:
        return "%2.2f " % num
    #
    def makeTermCollection(molecule_stoichiometries):
      """
      Formats a set of terms with stoichiometries.
      :param list-MoleculeStoichiometry:
      :return str:
      """
      term_collection = ''
      for m_s in molecule_stoichiometries:
        term = "%s%s" % (makeStoichiometryString(m_s), str(m_s.molecule))
        if len(term_collection) == 0:
          term_collection += term
        else:
          term_collection += " + " + term
      return term_collection
    #
    reactant_collection = makeTermCollection(self.reactants)
    product_collection = makeTermCollection(self.products)
    formula_str = ''
    reaction_str = "%s: %s -> %s" % (self.label,
        reactant_collection, product_collection)
    reaction_str = reaction_str + formula_str
    return reaction_str

  def reduceBySOMs(self):
    """
    Cancel reactants and products by SOMs
    SOMs is given by existing MESGraph.
    Return a new identifier. 
    :return str:
    """
    reactant_soms = {self.mesgraph.getNode(r.molecule) for r in self.reactants}
    product_soms = {self.mesgraph.getNode(p.molecule) for p in self.products}
    common_soms = list(reactant_soms.intersection(product_soms))
    if common_soms:
      for som in common_soms:
        reactants_in = collections.deque([ms for ms in
                                          self.reactants if
                                          self.mesgraph.getNode(ms.molecule)==som])
        reactants_out = [ms for ms in
                         self.reactants if
                         self.mesgraph.getNode(ms.molecule)!=som]
        products_in = collections.deque([ms for ms in
                                         self.products if
                                         self.mesgraph.getNode(ms.molecule)==som])
        products_out = [ms for ms in
                        self.products if
                        self.mesgraph.getNode(ms.molecule)!=som]
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
        if (len(self.reactants) > len(reactants)) | \
          (len(self.products) > len(products)):
          self.reactants = reactants
          self.products = products
      self.identifier = self.makeIdentifier()
    return self.identifier


class GAMESReport(object):

  def __init__(self, mesgraph, explain_threshold=20, errors=None):
    self.mesgraph = mesgraph
    self.errors = errors
    self.report_type_one_errors = NULL_STR
    self.explain_threshold = explain_threshold

  def getMoleculeEqualityPath(self, som, molecule1, molecule2):
    """
    Create an undirected graph between
    two molecules within a SOM
    and find the shortest path
    :param SOM som:
    :param Molecule mole1:
    :param Molecule mole2:
    :return PathComponents: som_path
    """   
    # molecule1 = mole1.name
    # molecule2 = mole2.name
    # construct undirected graph
    subg = nx.Graph()
    # here, every reaction is 1-1 reaction
    for reaction in list(som.reactions):
      node1 = reaction.reactants[0].molecule.name
      node2 = reaction.products[0].molecule.name
      if subg.has_edge(node1, node2):
        reaction_label = subg.get_edge_data(node1, node2)[cn.REACTION]
        # if reaction.label is not included in the attribute, add its label
        if reaction.label not in set(reaction_label):
          reaction_label = reaction_label + [reaction.label]
      else:
        reaction_label = [reaction.label]    
      subg.add_edge(node1, node2, reaction=reaction_label)
    path = [short_p for short_p in nx.shortest_path(subg, 
                                                    source=molecule1, 
                                                    target=molecule2
                                                    )
           ]
    som_path = []
    for idx in range(len(path)-1):
      edge_reactions = subg.get_edge_data(path[idx], path[idx+1])[cn.REACTION]
      som_path.append(cn.PathComponents(node1=path[idx], 
                                        node2=path[idx+1],
                                        reactions=edge_reactions
                                        )
                     )
    return som_path

  def getMoleculeEqualityPathReport(self, molecule_name1, molecule_name2, reaction_count, explain_details):
    """
    Print out shortest path between two molecules within a 
    same SOM. Molecule names (str) are arguments. 
    :param str molecule_name1:
    :param str molecule_name2:
    :param int reaction_count:
    :param bool explain_details:
    :return bool/str: path_report
    :return int: reaction_count
    """
    path_report = NULL_STR
    som1 = self.mesgraph.getNode(self.mesgraph.simple.getMolecule(molecule_name1))
    som2 = self.mesgraph.getNode(self.mesgraph.simple.getMolecule(molecule_name2))
    if som1 != som2:
      return False
    else:
      # in case molecule_name1 == molecule_name2, for example A -> A + B
      if molecule_name1 == molecule_name2:
      	if explain_details:
          path_report = path_report + "Clearly, %s %s %s\n" % (
              molecule_name1, cn.EQUAL, molecule_name2)
      else:
        som_path = self.getMoleculeEqualityPath(som1, molecule_name1, molecule_name2)
        for pat in som_path:
          if explain_details:
            path_report = path_report + "\n%s %s %s by reaction(s):\n" % (pat.node1, cn.EQUAL, pat.node2)
          for r in pat.reactions:
            som_reaction = self.mesgraph.simple.getReaction(r)
            reaction_count += 1
            path_report = path_report + "%d. %s\n" % (reaction_count, som_reaction.makeIdentifier(is_include_kinetics=False))
      return reaction_count, path_report

  def getMoleculeInequalityPathReport(self, molecule_name1, molecule_name2, reaction_names, reaction_count, explain_details):
  	"""
  	Print the reactions that infer inequality between molecules. 
  	Molecule name are given as arguments.
  	Especially, the mass of molecule 1 is less than that of molecule 2
  	:param str molecule_name1:
  	:param str molecule_name2:
  	:param list-str reaction_names:
  	:param int reaction_count:
  	:param bool explain_details
  	:return str: path_report
  	:return int: reaction_count
  	"""
  	path_report = NULL_STR
  	if explain_details:
  	  path_report = path_report + "%s %s %s by reaction(s):\n" % (molecule_name1, cn.LESSTHAN, molecule_name2)
  	for reaction_name in reaction_names:
  	  reaction_count += 1
  	  reaction = self.mesgraph.simple.getReaction(reaction_name)
  	  path_report = path_report + "%d. %s\n" % (reaction_count, reaction.makeIdentifier(is_include_kinetics=False))
  	return reaction_count, path_report

  def reportTypeOneError(self, type_one_errors, explain_details=False):
    """
    Generate report for Type I Errors. 
    Type I error occurs when there is a reaction
    that needs to add an arc between molecules,
    while the molecules are already included
    in the same SOM.
    :param list-PathComponents type_one_errors:
    :param bool explain_details:
    :return str: type_one_report
    :return list-int: error_num
    """
    report = NULL_STR
    error_num = []
    if len(type_one_errors) == 0:
      return report, error_num
    for pc in type_one_errors:
      mole1 = pc.node1
      mole2 = pc.node2
      reactions = pc.reactions
      reaction_count = 0
      reaction_count, equality_report = self.getMoleculeEqualityPathReport(mole1, mole2, reaction_count, explain_details)
      report = report + equality_report
      if explain_details:
        report = report + "\nHowever, "
      reaction_count, inequality_report = self.getMoleculeInequalityPathReport(mole1, mole2, reactions, reaction_count, explain_details)
      report = report + inequality_report
      report = report + "\n%s\n" % (PARAGRAPH_DIVIDER)
      error_num.append(reaction_count)
    report = report + "\n%s\n" % (REPORT_DIVIDER)
    return report, error_num

  def reportTypeTwoError(self, type_two_errors, explain_details=False):
  	"""
  	Generate report for Type II Errors.
  	Type II Error occurs when there is
  	a cycle between SOMs, which
  	should not happen.
  	:param list-(list-SOMs) type_two_errors:
  	:param bool explain_details:
  	:return str: type_two_report
  	:return list-int: error_num
  	"""
  	report = NULL_STR
  	error_num = []
  	if len(type_two_errors) == 0:
  	  return report, error_num
  	for cycle in type_two_errors:
  	  all_soms = set()
  	  reaction_count = 0
  	  report = report + "We detected a mass imbalance from the following reactions:\n"
  	  if explain_details:
  	    report = report + "\nThese uni-uni reactions created mass-equivalence.\n"
  	  for som in cycle:
  	  	all_soms.add(som)
  	  	for reaction in list(som.reactions):  	  	  	
  	  	  reaction_count += 1
  	  	  report = report + "\n%d. %s" % (reaction_count, reaction.makeIdentifier(is_include_kinetics=False))
  	  if explain_details:
  	  	report = report + "\n%s\n" % (PARAGRAPH_DIVIDER)
  	  	report = report + "The following reactions create mass-inequality.\n"
  	  inequality_reactions = []
  	  for som_pair in zip(cycle, cycle[1:] + [cycle[0]]):
  	    inequality_reactions = inequality_reactions + self.mesgraph.get_edge_data(som_pair[0], som_pair[1])[cn.REACTION]
  	  for r in inequality_reactions:
  	    reaction = self.mesgraph.simple.getReaction(r)
  	    reaction_count += 1
  	    report = report + "\n%d. %s" % (reaction_count, reaction.makeIdentifier(is_include_kinetics=False))
  	  if explain_details:
  	    report = report + "\n%s\n" % (PARAGRAPH_DIVIDER)
  	  if explain_details:
  	  	reaction_count = reaction_count - len(inequality_reactions) 	
  	  	# explain the mass equivalent pseudo reactions and the resulting sets
  	  	report = report + "Based on the reactions above, we have mass-equivalent pseudo reactions.\n"
  	  	for r in inequality_reactions:
  	  	  reaction = self.mesgraph.simple.getReaction(r)
  	  	  som_reaction = self.mesgraph.convertReactionToSOMReaction(reaction)
  	  	  reaction_count += 1
  	  	  report = report + "\n(pseudo %d.) %s" % (reaction_count, som_reaction.identifier)
  	  	report = report +  "\n%s\n" % (PARAGRAPH_DIVIDER)
  	  	report = report + "However, the above pseudo reactions imply the following inequalities:\n\n"
  	  	for som in cycle:
  	  	  report = report + "%s %s " % (som.identifier, cn.LESSTHAN)
  	  	report = report + "%s\n" % (cycle[0].identifier)
  	  	report = report + "\nThis indicates a mass conflict between reactions."
  	  error_num.append(reaction_count)
  	  report = report + "\n%s%s\n" % (PARAGRAPH_DIVIDER, PARAGRAPH_DIVIDER)
  	return report, error_num

  def convertOperationSeriesToReactionOperations(self, operation):
    """
    Convert an operation series to 
    a list of reaction operations.
    An 'operation' is either a row or column 
    of an operation matrix, 
    where both column and row indices are reactions. 
    :param pandas.Series operation:
    :return list-ReactionOperation: operations
    """
    operations = []
    # 
    nonzero_idx = np.array([idx for idx, val in enumerate(operation) if val != 0.0])
    nonzero_op = operation[nonzero_idx]
    for idx in range(len(nonzero_op)):
      reaction_op = ReactionOperation(reaction=nonzero_op.index[idx],
      	                              operation=nonzero_op[idx]
      	                              )
      operations.append(reaction_op)
    return operations

  def getOperationMatrix(self):
  	"""
  	Return an operation matrix (pandas DataFrame)
  	on the transposed stoichiometry matrix
  	:return pandas.DataFrame: operation_df
  	"""
  	operation_df = None
  	if self.mesgraph.lower_inverse is None:
  	  pass
  	elif self.mesgraph.rref_operation is None:
  	  operation_df = self.mesgraph.lower_inverse
  	else:
  	  operation_df = self.mesgraph.rref_operation.dot(self.mesgraph.lower_inverse)
  	return operation_df

  def getResultingSeries(self, reaction_label):
    """
    Return a reaction series, that is, a column
    from a reduced stoichiometry matrix.
    :param str reaction_label:
    :return False/pandas.Series: result_series
    """
    if type(reaction_label) != str:
      return False
    else: 
      if self.mesgraph.rref_df is None:
        result_series = self.mesgraph.echelon_df[reaction_label]
      else:
        result_series = self.mesgraph.rref_df[reaction_label]
      return result_series

  def getOperationStoichiometryMatrix(self, reaction_operations):
  	"""
  	Create a stoichiometry matrix of reactions 
  	related to operations.
  	:param list-ReactionOperation reaction_operations:
  	:return pandas.DataFrame: stoichiometry_df
  	"""
  	reactions = []
  	species = set()
  	for op in reaction_operations:
  	  reaction = self.mesgraph.simple.getReaction(op.reaction)
  	  reactions.append(reaction)
  	  species = species.union({r.molecule.name for r in reaction.reactants})
  	  species = species.union({p.molecule.name for p in reaction.products})
  	operations = np.array([val.operation for val in reaction_operations])
  	stoichiometry_df = self.mesgraph.getStoichiometryMatrix(reactions, list(species))
  	return stoichiometry_df

  def getInferredReaction(self, reaction_operations):
  	"""
  	Create an inferred reaction from reaction_operations.
  	An inferred reaction is a linear combination of reactions.
  	:param list-ReactionOperation reaction_operations:
  	:return SimplifiedReaction: inferred_reaction
  	"""
  	INFERRED_REACTION = "Inferred Reaction"
  	stoichiometry_df = self.getOperationStoichiometryMatrix(reaction_operations)
  	reaction_index = [op.reaction for op in reaction_operations]
  	operation_series = pd.Series([val.operation for val in reaction_operations], index=reaction_index)
  	resulting_reaction = stoichiometry_df.dot(operation_series)
  	# Create a simplified reaction based on the resulting_reaction
  	reactants = []
  	products = []
  	for ms in resulting_reaction.iteritems():
  	  if abs(ms[1]) < TOLERANCE:
  	  	continue
  	  elif ms[1] > 0:
  	  	products.append(MoleculeStoichiometry(molecule=self.mesgraph.simple.getMolecule(ms[0]),
  	  	                                      stoichiometry=ms[1]))
  	  else:
  	  	reactants.append(MoleculeStoichiometry(molecule=self.mesgraph.simple.getMolecule(ms[0]),
  	  	                                      stoichiometry=abs(ms[1])))
  	inferred_reaction = SimplifiedReaction(reactants,
  		                                   products,
  		                                   NULL_STR,
  		                                   self.mesgraph)
  	inferred_reaction.reduceBySOMs()
  	return inferred_reaction

  def reportReactionsInSOM(self, som, reaction_count=0):
    """
    Generate report on reactions created a SOM,
    in order to demonstrate how the SOM was constructed.
    :param som SOM:
    :param int reaction_count:
    :return str: report
    :return int: reaction_count
    """
    report = cn.NULL_STR
    if not som.reactions:
      return report, reaction_count
    reactions = list(som.reactions)
    for r in reactions:
      reaction_count += 1
      report = report + "%d. %s\n" % (reaction_count,
                                      r.makeIdentifier(is_include_kinetics=False))
    return report, reaction_count

  def reportTypeThreeError(self, type_three_errors, explain_details=False):
    """
    Generate a report for Type III errors.
    A Type III error occurs when there is 
    a 1-1 SOMReaction, while there is a preexisting
    MESGraph node between them. 
    :param list-SOMReaction type_three_errors:
    :param bool explain_details:
    :return False/str: report
    :return list-int: error_num
    """
    report = NULL_STR
    error_num = []
    if len(type_three_errors) == 0:
      return report, error_num
    for type3_error in type_three_errors:
      reaction_count = 0
      if type3_error.category != cn.REACTION_1_1:
        print("This canot be a type three error!")
        report = False
        break
      else:
        reaction_label = type3_error.label
        reactant_som = type3_error.reactants[0].som
        product_som = type3_error.products[0].som
        operation_df = self.getOperationMatrix()
        operation_series = operation_df.T[reaction_label]
        reaction_operations = self.convertOperationSeriesToReactionOperations(operation_series)
        # if the number of elements exceeds the threshold, do not explain details
        if len(reaction_operations) > self.explain_threshold:
          explain_details = False
        inferred_reaction = self.getInferredReaction(reaction_operations)
        inferred_som_reaction = self.mesgraph.convertReactionToSOMReaction(inferred_reaction)
        if inferred_som_reaction.getCategory() != cn.REACTION_1_1:
          explain_details = False
        reactant_som = inferred_som_reaction.reactants[0].som
        product_som = inferred_som_reaction.products[0].som
        inequality_reactions = []
        if self.mesgraph.has_edge(reactant_som, product_som):
          inequality_reactions = inequality_reactions + self.mesgraph.get_edge_data(reactant_som, product_som)[cn.REACTION]
        if self.mesgraph.has_edge(product_som, reactant_som):
          inequality_reactions = inequality_reactions + self.mesgraph.get_edge_data(product_som, reactant_som)[cn.REACTION]
        #
        report = report + "We detected a mass imbalance from the following reactions:\n"
        soms = set()
        som_reactions = []
        for r in reaction_operations:
          reaction = self.mesgraph.simple.getReaction(r.reaction)
          reaction_count += 1
          report = report + "\n%d. %s" % (reaction_count, reaction.makeIdentifier(is_include_kinetics=False))
          som_reaction = self.mesgraph.convertReactionToSOMReaction(reaction)
          som_reactions.append(som_reaction)
          soms = soms.union({r.som for r in som_reaction.reactants})
          soms = soms.union({p.som for p in som_reaction.products})
        # calculated the reported number of reactions that constitute a type III error
        error_num.append(reaction_count + len(inequality_reactions))
        if explain_details:
          report = report + "\n\n%s%s\n" % ("-"*NUM_STAR, PARAGRAPH_DIVIDER)
          report = report + "These uni-uni reactions created mass-equivalence.\n" 
          report = report + "(The chemical species within a curly bracket have the same atomic mass.)\n"
        else:
          report = report + "\n"  	  	
        for som in soms:
          if explain_details and som.reactions:
            report = report + "\n%s is inferred by:\n" % som.makeId()
          sub_report, reaction_count = self.reportReactionsInSOM(som, reaction_count)
          report = report + sub_report
        if explain_details:
          report = report + "%s\n" % (PARAGRAPH_DIVIDER)
          report = report + "These multi-uni reactions created mass-inequality.\n\n"  	  
        for r in inequality_reactions:
          reaction = self.mesgraph.simple.getReaction(r)
          reaction_count += 1
          report = report + "%d. %s\n" % (reaction_count, reaction.makeIdentifier(is_include_kinetics=False))
        reaction_count = reaction_count - len(inequality_reactions)
        if explain_details:
          report = report + "%s\n" % (PARAGRAPH_DIVIDER)
          report = report + "Based on the reactions above, we have mass-equivalent pseudo reactions.\n"
          pseudo_reaction_count = 0
          for sr in som_reactions:
            pseudo_reaction_count += 1
            report = report + "\n(pseudo %d.) %s" % (pseudo_reaction_count, sr.identifier)
          report = report +  "\n%s\n" % (PARAGRAPH_DIVIDER)
          report = report + "An operation between pseudo reactions:\n"
          report = report + "\n%.2f * %s" % (reaction_operations[0].operation, reaction_operations[0].reaction)
          for ro in reaction_operations[1:]:
            if ro.operation < 0:
              report = report + " - "
            else:
              report = report + " + "
            report = report + "%.2f * %s\n" % (abs(ro.operation), ro.reaction)
          report = report + "\n\nwill result in a uni-uni reaction:\n"
          report = report + "\n%s\n" % (inferred_som_reaction.identifier)
          report = report + "\n\nmeaning %s and %s have equal mass.\n" % (reactant_som, product_som)
          report = report +  "%s\n" % (PARAGRAPH_DIVIDER)
          report = report + "However, the following mass-equivalent pseudo reaction(s):\n"
          for r in inequality_reactions:
            reaction = self.mesgraph.simple.getReaction(r)
            reaction_count += 1
            report = report + "\n%d. %s" % (reaction_count, reaction.makeIdentifier(is_include_kinetics=False))
            som_reaction = self.mesgraph.convertReactionToSOMReaction(reaction)
            report = report + "\n(pseudo %d.) %s" % (reaction_count, som_reaction.identifier)
          report = report + "\n\nincidates the masses of %s and %s are unequal.\n" % (reactant_som, product_som)
          #
          report = report +  "\n%s%s\n" % (PARAGRAPH_DIVIDER, PARAGRAPH_DIVIDER)
        report = report + "\n%s\n" % (REPORT_DIVIDER) 
    return report, error_num

  def reportEchelonError(self, echelon_errors, explain_details=False):
    """
    Generate a report for echelon errors, i.e.
    mass balance errors from LU decomposition/RREF
    and store in self.report_echelon_errors.
    The operation_df is either mesgraph.lower_inverse
    if LU decomposition yielded echelon_errors,
    and rref_df.dot(lower_inverse) if RREF created errors. 
    :param list-SOMReaction echelon_errors:
    :param bool explain_details:
    :return str: report
    :return list-int: error_num
    """
    report = NULL_STR
    error_num = []
    if len(echelon_errors) == 0:
      return report, error_num
    operation_df = self.getOperationMatrix()
    error_report = NULL_STR
    for reaction in echelon_errors:
      reaction_count = 0
      reaction_label = reaction.label
      operation_series = operation_df.T[reaction_label]
      result_series = self.getResultingSeries(reaction_label)
      reaction_operations = self.convertOperationSeriesToReactionOperations(operation_series)
      # if the number of elements exceeds the threshold, do not explain details
      if len(reaction_operations) > self.explain_threshold:
        explain_details = False
      inferred_reaction = self.getInferredReaction(reaction_operations)
      inferred_som_reaction = self.mesgraph.convertReactionToSOMReaction(inferred_reaction)
      #
      report = report + "\nWe detected a mass imbalance\n%s\n" % inferred_reaction.identifier
      report = report + "\nfrom the following isolation set.\n"
      #
      nonzero_idx = np.array([idx for idx, val in enumerate(result_series) if val != 0.0])
      nonzero_result_series = result_series[nonzero_idx]
      #
      # part 1: reactions that caused mass balance errors
      reaction_operations = self.convertOperationSeriesToReactionOperations(operation_series)
      reported_reactions = [r.reaction for r in reaction_operations]
      reported_som_reactions = []
      for r in reported_reactions:
      	reaction = self.mesgraph.simple.getReaction(r)
      	reaction_count += 1
      	report = report + "\n%d. %s" % (reaction_count, reaction.makeIdentifier(is_include_kinetics=False))
      error_num.append(reaction_count)
      one_side = "--undetermined--"
      if inferred_som_reaction.reactants==[]:
        one_side = "reactant"
      elif inferred_som_reaction.products==[]:
        one_side = "product"
      if one_side == "--undetermined--":
        explain_details = False
      #
      # part 2: SOMs that were canceled by the operation
      canceled_soms = set()
      for r in reported_reactions:
      	sr = self.mesgraph.convertReactionToSOMReaction(self.mesgraph.simple.getReaction(r))
      	reported_som_reactions.append(sr)
      	canceled_soms = canceled_soms.union({r.som for r in sr.reactants})
      	canceled_soms = canceled_soms.union({p.som for p in sr.products})
      canceled_soms = canceled_soms.difference({r.som for r in inferred_som_reaction.reactants})
      canceled_soms = canceled_soms.difference({p.som for p in inferred_som_reaction.products})
      #
      if explain_details:
      	report = report + "\n\n%s%s\n" % ("-"*NUM_STAR, PARAGRAPH_DIVIDER)
      	report = report + "These uni-uni reactions created mass-equivalence.\n"
      	report = report + "(The chemical species within a curly bracket have the same atomic mass.)\n"
      else:
        report = report + "\n"
      for som in canceled_soms:
        if explain_details and som.reactions:
          report = report + "\n%s is inferred by:\n" % som.makeId()
        sub_report, reaction_count = self.reportReactionsInSOM(som, reaction_count)
        report = report + sub_report
      if explain_details:
        report = report + "%s\n" % (PARAGRAPH_DIVIDER)
        report = report + "Based on the uni-uni reactions above, we create mass-equivalent pseudo reactions.\n"
        pseudo_reaction_count = 0
        for sr in reported_som_reactions:
          pseudo_reaction_count += 1
          report = report + "\n(pseudo %d.) %s" % (pseudo_reaction_count, sr.identifier)
        report = report +  "\n%s\n" % (PARAGRAPH_DIVIDER)
        report = report + "An operation between the pseudo reactions:\n"
        report = report + "%.2f * %s" % (reaction_operations[0].operation, reaction_operations[0].reaction)
        for ro in reaction_operations[1:]:
          if ro.operation < 0:
            report = report + " - "
          else:
            report = report + " + "
          report = report + "%.2f * %s" % (abs(ro.operation), ro.reaction)
        # one_side is determined in lines 589-593
        report = report + "\n\nwill result in empty %s with zero mass:\n" % (one_side)
        report = report + "\n%s\n" % (inferred_som_reaction.identifier)
        report = report +  "\n%s%s\n" % (PARAGRAPH_DIVIDER, PARAGRAPH_DIVIDER)
      report = report + "\n%s\n" % (REPORT_DIVIDER)
    return report, error_num

  def reportCancelingError(self, canceling_errors, explain_details=False):
    """
    Generate a report for canceling errors.
    A canceling error occurs when a som_reaction has
    an imbalanced net stoichiometry.
    :param list-SOMReaction:
    :return str: report
    :return list-int: error_num
    """
    report = NULL_STR
    error_num = []
    if len(canceling_errors) == 0:
      return report, error_num
    for error in canceling_errors:
      reaction_count = 0
      label = error.label
      reaction = self.mesgraph.simple.getReaction(label)
      simplified_reaction = SimplifiedReaction(reaction.reactants, reaction.products, label, self.mesgraph)
      simplified_reaction.reduceBySOMs()
      som_reaction = self.mesgraph.convertReactionToSOMReaction(reaction)
      som_reactants = {r.som for r in som_reaction.reactants}
      som_products = {p.som for p in som_reaction.products}
      canceled_soms = list(som_reactants.intersection(som_products))
      reaction_count += 1
      report = report + "\n%d. %s" % (reaction_count, reaction.makeIdentifier(is_include_kinetics=False))
      for som in canceled_soms:
        if len(som.molecules) == 1 and explain_details:
          molecule = list(som.molecules)[0]
          report = report + "\n*%s is a common chemical species in reactants and products, so can be canceled\n" % (molecule.name)
        for reaction in list(som.reactions):
          reaction_count += 1
          reactant = reaction.reactants[0].molecule
          product = reaction.products[0].molecule
          report = report + "\n%d. %s" % (reaction_count, reaction.makeIdentifier(is_include_kinetics=False))
          if explain_details:
            report = report + "\n*%s and %s have the same mass according to the above reaction\n" % (reactant.name, product.name)
      if error.reactants==[]:
        one_side = "reactant"
      elif error.products==[]:
        one_side = "product"
      report = "We detected a mass imbalance\n: %s\n\nfrom the following isolation set:\n" % (simplified_reaction.identifier) + report 
      error_num.append(reaction_count)      
      report = report + "\n%s%s\n" % (PARAGRAPH_DIVIDER, PARAGRAPH_DIVIDER)
    report = report + "\n%s\n" % (REPORT_DIVIDER)
    return report, error_num








