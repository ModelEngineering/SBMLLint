"""
Sample program to illustrate how to extract reactants, products, and kinetics expressions.
1. Create objects
2. Handle stoichiometry
"""
import sys
import os.path
import libsbml

def getModel(filename):
  """
  :param str filename: input file with the xml model
  :return libsbml.Model"
  """
  reader = libsbml.SBMLReader()
  document = reader.readSBML(filename)
  if (document.getNumErrors() > 0):
      print("Encountered the following SBML errors:" + "\n");
      document.printErrors();
      exit(-1)
  model = document.getModel()
  import pdb; pdb.set_trace()
  return model

def getReactions(model):
  """
  :param libsbml.Model:
  :return list-of-reactions
  """
  num = model.getNumReactions()
  return [model.getReaction(n) for n in range(num)]

def getParameters(model):
  """
  :param libsbml.Model:
  :return list-of-reactions
  """
  return [model.getParameter(n) for n in range(model.getNumParameters())]

def getReactants(reaction):
  """
  :param libsbml.Reaction:
  :return list-of-libsbml.SpeciesReference:
  """
  return [reaction.getReactant(n) for n in range(reaction.getNumReactants())]

def getProducts(reaction):
  """
  :param libsbml.Reaction:
  :return list-of-libsbml.SpeciesReference:
  """
  return [reaction.getProduct(n) for n in range(reaction.getNumProducts())]

def getReactionString(reaction):
  reaction_str = ''
  base_length = len(reaction_str)
  for reference in getReactants(reaction):
    if len(reaction_str) > base_length:
      reaction_str += " + " + reference.species
    else:
      reaction_str += reference.species
  reaction_str += "-> "
  base_length = len(reaction_str)
  for reference in getProducts(reaction):
    if len(reaction_str) > base_length:
      reaction_str += " + " + reference.species
    else:
      reaction_str += reference.species
  kinetics_terms = getReactionKineticsTerms(reaction)
  reaction_str += "; " + ", ".join(kinetics_terms)
  return reaction_str

def getReactionKineticsTerms(reaction):
  """
  Gets the terms used in the kinetics law for the reaction
  :param libsbml.Reaction
  :return list-of-str: names of the terms
  """
  terms = []
  law = reaction.getKineticLaw()
  if law is not None:
    math = law.getMath()
    asts = [math]
    while len(asts) > 0:
      this_ast = asts.pop()
      if this_ast.isName():
        terms.append(this_ast.getName())
      else:
        pass
      num = this_ast.getNumChildren()
      for idx in range(num):
        asts.append(this_ast.getChild(idx))
  return terms
    
def getKineticsParametersAndSpecies(reaction):
  """
  :param Reaction reaction:
  :return list-of-str, list-of-str:
  """

  

def main(args):
  if (len(args) != 2):
      print("\n" + "Usage: printMath filename" + "\n" + "\n");
      return 1;
  filename = args[1];
  reader = libsbml.SBMLReader()
  document = reader.readSBML(filename)
  if (document.getNumErrors() > 0):
      print("Encountered the following SBML errors:" + "\n");
      document.printErrors();
      exit(-1)
  model = document.getModel()
  reactions = getReactions(model)
  parameters = getParameters(model)
  for reaction in reactions:
    print getReactionString(reaction)

if __name__ == '__main__':
  main(sys.argv)  
