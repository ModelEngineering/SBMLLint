# -*- coding: utf-8 -*-
"""Doctrings on StoichiometryMatrix.py.

The StoichiometryMatrix class creates a stoichiometry matrix class of a model and tests
whether it is consistent. The columns of stoichiometry matrix represent reactions
and the rows represent metabolites. 

A stoichiometry matrix is 'consistent' if, for a matrix N and a metabolite mass vector m,
we can find positive vector m such that N.T*m = 0. It is inconsistent if we cannot find
such a solution. The isConsistent method checks whether the stoichiometry matrix is 
consistent or not. 

Example:
    You should already have an SBML model loaded, which we call 'model'.

    >> import StoichiometryMatrix as ss
    >> sbml_mat = ss.StoichiometryMatrix(model)
    >> sbml_mat.isConsistent()

Notes:
    * If Tellurium is updated with a new stoichiometry matrix in future, we may simply use it instead

"""

# Help from https://gist.github.com/lukauskas/d1e30bdccc5b801d341d


import constants as ct
import numpy as np
import pandas as pd
from scipy.optimize import linprog
import tesbml

class StoichiometryMatrix:
    """Creates a full stoichiometry matrix from a model.

    This class takes one input argument, an SBML model, and creates
    a MassFlowGraph and a FullGraph. The user can check the model's
    stoichiometric consistency using isStoichiometryConsistent(). 

    Attributes:
        model (SBML model): An SBML model. The argument needed to create a ModelGraph object. 
        connected (str): Connectivity status of a model.

        stoichiometry_matrix: A full stoichiometry matrix. By convention, rows represent species
        and columns represent reactions. 

    """    
    def __init__(self, model):
        self.model = model

        if self.model == None:
            raise TypeError("Model doesn't exist")
    
    def buildMatrix(self):
        """Creates a full stoichiometry matrix from a model.
        
        We assume each reaction has at least one reactant and one product.
        if either is undefined as EmptySet, we create a new boundary species 
        with stoichiometry +/-1.0. This full stoichiometry matrix is used to
        determine stoichiometric inconsistency using linear programming.
        
        Args:
            No specific arguments. 

        Returns:
            Stoichiometry matrix (pandas Dataframe)    
        """
        
        list_reactions = [reaction.getId() for reaction in self.model.getListOfReactions()]
        set_reactions = {reaction.getId() for reaction in self.model.getListOfReactions()}
        set_species = {species.getId() for species in self.model.getListOfSpecies()}.difference({ct.EMPTYSET})
        list_species = list(set_species)

        stoichiometry_matrix = {}
        num_reactions = self.model.getNumReactions()

        # Rows: species, Columns: reactions
        for reaction_idx in range(num_reactions):

            reaction = self.model.getReaction(list_reactions[reaction_idx])
            list_reactants = [reactant.getSpecies() for reactant in reaction.getListOfReactants()]
            list_products = [product.getSpecies() for product in reaction.getListOfProducts()]

            reactants = {reactant.getSpecies(): reactant.getStoichiometry() for reactant in reaction.getListOfReactants()}
            if (reactants == {}) | (reactants.keys() == {ct.EMPTYSET}):
                BDRY_REACTANT_NAME = reaction.getId()+'_bdry_rct'
                reactants = {BDRY_REACTANT_NAME: 1.0}
                list_species.append(BDRY_REACTANT_NAME)

            products = {product.getSpecies(): product.getStoichiometry() for product in reaction.getListOfProducts()}
            if (products == {}) | (products.keys() == {ct.EMPTYSET}):
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

        self.stoichiometry_matrix = full_stoichiometry_matrix   
        return full_stoichiometry_matrix

    
    def isConsistent(self):
        """ Runs linear programmming to determine inconsistency. 
        
        Status 0 means optimization terminated successfully, 
        and 2 means the value was infeasible, i.e. the model is
        inconsistent.

        Args:
            No specific arguments. 

        Returns:
            True if consistent, False if not consistent
        """

        s_matrix_t = self.stoichiometry_matrix.T

        # number of metabolites
        nmet = s_matrix_t.shape[0]
        # number of reactions
        nreac = s_matrix_t.shape[1]

        b = np.zeros(nmet)
        c = np.ones(nreac)

        # linear programming. c is constraint (here, zero), 
        # and b is possible values for metabolite vector. 
        res = linprog(c, A_eq=s_matrix_t, b_eq=b, bounds=(1, None))

        if res.status == 0:
            self.consistent = True
        else:
            self.consistent = False

        return self.consistent
















