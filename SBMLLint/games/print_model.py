# -*- coding: utf-8 -*-
"""Doctrings on print_model.py.

This script has a simple function named print_model,
which prints out all reactions in a given model.
The format is as follows:
'reaction' 'reaction index': 'reactants' '->' 'products'

The goal is to help examine models, for example
to manually verify whether a model is stoichiometrically inconsistent.

"""



def print_model(model):
    """Print reactants and products of all reactions in a model.


    Args:
        model (sbml_model): An sbml model to print out

    Returns:
        None: Just prints out the reactions
        
    """

    # print model name
    print(model)

    # run for each reaction
    for reaction_idx in range(model.getNumReactions()):
        reaction = model.getReaction(reaction_idx)

        #stg = "reaction " + str(reaction_idx+1) + ": "
        stg = reaction.getId() + ": "

        reactants_list = [reactant.getSpecies() for reactant in reaction.getListOfReactants()]
        products_list = [product.getSpecies() for product in reaction.getListOfProducts()]
        stg = stg + ' + '.join(reactants_list) + " -> " + ' + '.join(products_list) + ';'
        
        print(stg)
