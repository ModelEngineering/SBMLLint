<img src="https://travis-ci.org/ModelEngineering/SBMLLint.svg?branch=master" width="100"/>

# SBMLlint

## Problem Addressed

Many biological models are based on chemical reactions. For example, glycolysis, arguably the most widely exercised metabolic pathway in biology, begins by transforming the reactants glucosue (Glu) and adenosine triphosphate (ATP) into the products glucose 6-phosphate (GluP) and adenosine diphosphate (ADP), or: ``Glu + ATP -> GluP + ADP``. Examples of biological modeling techniques that rely on reactions include flux balance analysis and kinetics models.

Today's biological models typically consist of tens to hundreds of reactions. Even with this modest level of complexity, it is easy to make mistakes. For example, a **mass balance error** occurs if the total mass of the reactants differs from the total mass of the products. With the advent of high throughput laboratory techniques, the complexity of models will grow rapidly. As a point of comparison, the complexity of software has grown from hundereds of lines of code in the 1960s to tens of millions of lines of code for software such as linux and the apache web server.

Because of this huge growth in the complexity, software engineers developed sophisticated tools to detect errors in codes *statically*, before any statement is executed. For example, the ``pylint`` tool analyzes ``python`` source codes to determine if a variable is referenced before a value is assigned to it. The term **linter** is used for a tool that does static analysis of source codes.

## The Tool

``SBMLLint`` is a tool that lints reactions. The initial focus is detecting mass balance errors. The tool takes as input a model expressed in either SBML (Systems Biology Markup Language) or the antimony language.

``SBMLLint`` implements two algorithms for linting reactions. The first, ``structured_names``, requires the modeller to give the tool hints by naming molecules in terms of their underlying moieties (sub-parts). For example, ``ATP`` would be written as ``A_P_P_P`` to indicate that there is one adenosine and three phosphates. And ``GluP`` would be written as ``Glu_P`` to indicate that it has one phosphate. Thus, the above reaction is written as ``Glu + A_P_P_P -> Glu_P + A_P_P``. ``structured_names`` checks that the count of each moiety in the reactants is the same as the count of each moiety in the products. Although ``structured_names`` places a burden on the modeller, we note that about 20% of the models in the BioModels repository already use names structured in the manner required by this tool. 

The second algorithm, ``games``, does not impose any requirements on the structure of the molecule names. However, ``games`` checks for a weaker condition called *sthoichiometry inconsistency*. A collections of reactions is stoichiometrically inconsistent if the set of reactions infers that a molecule has a mass less than it's own mass. To illustrate this, consider the reactions ``A -> B + C`` and ``C -> A``. The first reaction implies that the mass of ``A`` is greater than the mass of ``C``. But the second reaction implies that ``A`` and ``C`` have the same mass.

## Installation
``SBMLLint`` is currently under development and so is not availabe for installation. An installable version is expected by April, 2019.

## Future Directions

Features under consideration for future versions include:

- Use of parameters or species that are uninitialized (and
not present on the right hand side of a reaction)
- Correctly constructed kinetics laws
- Parameters that are unreferenced
- Mass conservation

