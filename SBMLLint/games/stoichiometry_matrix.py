"""Stoichiometry Matrix."""

from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule, MoleculeStoichiometry
from SBMLLint.common.reaction import Reaction
from SBMLLint.games.som import SOM
from SBMLLint.common.simple_sbml import SimpleSBML

import itertools
import numpy as np
import pandas as pd
from scipy.optimize import linprog

class StoichiometryMatrix(object):
  """
  Creates a full stoichiometry matrix from simpleSBML 
  that correctly incorporates boundary species 
  and use linear programming to determine stoichiometric consistency.
  """    
  def __init__(self, simple=None):
    self.reactions = self._getReactions(simple)
    self.matrix = self._getStoichiometryMatrix()
    self.consistent = None

  def _getReactions(simple):
    """
    Choose non-boundary reacetions
    :param SimpleSBML simple:
    :return list-Reaction:
    """
    reactions = []
    for reaction in simple.reactions:
      if reaction.category != cn.REACTION_BOUNDARY:
        reacitons.append(reaction)
    return reactions

  def _getStoichiometryMatrix(self):
    """
    Creates a full stoichiometry matrix
    using non-boundary reactions.
    :return pd.DataFrame:
    """
        
        list_reactions = [reaction.getId() for reaction in self.model.getListOfReactions()]
        set_reactions = {reaction.getId() for reaction in self.model.getListOfReactions()}
        set_species = {species.getId() for species in self.model.getListOfSpecies()}.difference({cn.EMPTYSET})
        list_species = list(set_species)

        stoichiometry_matrix = {}
        num_reactions = self.model.getNumReactions()

        # Rows: species, Columns: reactions
        for reaction_idx in range(num_reactions):

            reaction = self.model.getReaction(list_reactions[reaction_idx])
            list_reactants = [reactant.getSpecies() for reactant in reaction.getListOfReactants()]
            list_products = [product.getSpecies() for product in reaction.getListOfProducts()]
            reactants = {reactant.getSpecies(): reactant.getStoichiometry() for reactant in reaction.getListOfReactants()}
            if (reactants == {}) | (reactants.keys() == {cn.EMPTYSET}):
                BDRY_REACTANT_NAME = reaction.getId()+'_bdry_rct'
                reactants = {BDRY_REACTANT_NAME: 1.0}
                list_species.append(BDRY_REACTANT_NAME)

            products = {product.getSpecies(): product.getStoichiometry() for product in reaction.getListOfProducts()}
            if (products == {}) | (products.keys() == {cn.EMPTYSET}):
                BDRY_PRODUCT_NAME = reaction.getId()+'_bdry_pdt'
                products = {BDRY_PRODUCT_NAME: 1.0}
                list_species.append(BDRY_PRODUCT_NAME)

            for species_index, species_node in enumerate(list_species):
                net_stoichiometry =  int(products.get(species_node, 0)) - int(reactants.get(species_node, 0))
                stoichiometry_matrix[species_index, reaction_idx] = net_stoichiometry

        nonzero_stoichiometry_matrix = {key:itm for key,itm in stoichiometry_matrix.items() if itm!=0}

        full_stoichiometry_matrix = pd.DataFrame(0, index=list_species, columns=list_reactions)
        for key in nonzero_stoichiometry_matrix.keys():
            full_stoichiometry_matrix.iloc[key] = nonzero_stoichiometry_matrix.get(key)
        
        # delete remaining nonzero rows; i.e. species in getListOfSpecies but not in actual reactions 
        result_matrix = full_stoichiometry_matrix.loc[(full_stoichiometry_matrix!=0).any(axis=1)]

        self.stoichiometry_matrix = result_matrix   
        return full_stoichiometry_matrix

    # def isConsistent(self):
    #     """ Runs linear programmming to determine inconsistency. 
        
    #     Status 0 means optimization terminated successfully, 
    #     and 2 means the value was infeasible, i.e. the model is
    #     inconsistent.

    #     Args:
    #         None. 

    #     Returns:
    #         True if consistent, False if not consistent
    #     """
    #     s_matrix_t = self.stoichiometry_matrix.T
    #     # number of metabolites
    #     nmet = s_matrix_t.shape[0]
    #     # number of reactions
    #     nreac = s_matrix_t.shape[1]
    #     #
    #     b = np.zeros(nmet)
    #     c = np.ones(nreac)
    #     # linear programming. c is constraint (here, zero), 
    #     # and b is possible values for metabolite vector. 
    #     res = linprog(c, A_eq=s_matrix_t, b_eq=b, bounds=(1, None))
    #     if res.status == 0:
    #         self.consistent = True
    #     else:
    #         self.consistent = False
    #     #
    #     return self.consistent
















