"""Set of Molecules(SOM)."""

from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common.reaction import Reaction
from SBMLLint.common.simple_sbml import SimpleSBML

BRACKET_OPEN = "{"
BRACKET_CLOSE = "}"

class SOM(object):
    soms = []  # All SOMs. 
    def __init__(self, molecule):
        """
        :param molecules Molecule instances:
        :param reactions Reaction instances:
        """
        self.molecules = {molecule}
        self.reactions = set()
        self.identifier = self.makeId()
        self.__class__.addSOM(self)

    def __repr__(self):
        return self.identifier        
        
    def makeId(self):
        """
        Creates an identifier for the SOM to uniquely
        identifies its elements.
        :return str:
        """
        def joinMoleculeNames(molecules):
          names = [m.name for m in molecules]
          names.sort()
          return ', '.join(names)
        #
        identifier = "%s%s%s" % (
            BRACKET_OPEN, 
            joinMoleculeNames(list(self.molecules)),
            BRACKET_CLOSE
            )
        return identifier
        
    @classmethod    
    def addSOM(cls, new_som):
        if any([new_som.molecules.intersection(s.molecules) for s in cls.soms]):
          pass
        else:
          cls.soms.append(new_som)
    
    @classmethod
    def findSOM(cls, molecule):
        """
        Find the SOM that contains molecule
        and returns SOM
        :param Molecule molecule
        :return SOM
        """    
        for som in cls.soms:
            for m in som.molecules:
                if molecule.name == m.name:
                    return som
            
    @classmethod
    def merge(cls, reaction):
        """
        Merges two SOMs using a UniUni reaction 
        and updates cls.soms
        :param Reaction reaction
        :return SOM merged som
        """ 
        if reaction.category != cn.REACTION_1_1:
            raise AttributeError("This reaction cannot merge. You need a 1-1 reacton")
            pass
        som1 = cls.findSOM(reaction.reactants[0].molecule)
        som2 = cls.findSOM(reaction.products[0].molecule)
        
        if som1 == som2:
            pass
        else: 
            som1.molecules = som1.molecules.union(som2.molecules)
            som1.reactions.add(reaction)
            som1.identifier = som1.makeId()
            cls.soms.remove(som2)
            return som1
        
    @classmethod
    def reduce(cls, reaction):
        """
        Reduces reaction using existing cls.soms
        param Reaction reaction
        return reduced reaction if reduced
        return False if reaction is not reducible
        """        
        # flag that will show whether the reaction was reduced
        reduced = False
        # Quit if reaction is not MultiMulti
        if reaction.category != cn.REACTION_n_n:
            return False
        
        def getIndex(ms_list, som):
            for key, ms in enumerate(ms_list):
                if (ms.molecule in som.molecules) & (ms.stoichiometry != 0):
                    return key
            return len(ms_list)        

        def getNamedTuple(tup, stoichiometry):
            return cn.MoleculeStoichiometry(molecule = tup.molecule,
                                           stoichiometry = stoichiometry)        
        
        for som in cls.soms:
            
            rct_index = getIndex(reaction.reactants, som)
            pdt_index = getIndex(reaction.products, som)

            while (rct_index<len(reaction.reactants)) & (pdt_index<len(reaction.products)):

                reactant = reaction.reactants[rct_index]
                product = reaction.products[pdt_index]

                if reactant.stoichiometry > product.stoichiometry:
                    rct_stoichiometry = reactant.stoichiometry - product.stoichiometry
                    pdt_stoichiometry = 0.0

                elif reactant.stoichiometry < product.stoichiometry:
                    rct_stoichiometry = 0.0
                    pdt_stoichiometry = product.stoichiometry - reactant.stoichiometry

                else:
                    rct_stoichiometry = 0.0
                    pdt_stoichiometry = 0.0

                reaction.reactants[rct_index] = getNamedTuple(reaction.reactants[rct_index], \
                                                              rct_stoichiometry)
                reaction.products[pdt_index] = getNamedTuple(reaction.products[pdt_index], \
                                                             pdt_stoichiometry) 
                reduced = True

                rct_index = getIndex(reaction.reactants, som)
                pdt_index = getIndex(reaction.products, som)
        
        if reduced: 
            # Update reaction 
            for ms in reaction.reactants:
                if ms.stoichiometry == 0.0:
                    reaction.reactants.remove(ms)
            for ms in reaction.products:
                if ms.stoichiometry == 0.0:
                    reaction.products.remove(ms)

            reaction.identifier = reaction.makeId()
            reaction.category = reaction._getCategory()     
            return reaction
        else:
            return reduced
