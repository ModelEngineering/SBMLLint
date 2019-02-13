<img src="https://travis-ci.org/ModelEngineering/SBMLLint.svg?branch=master" width="100"/>

# SBMLlint
Model checker for [SBML](http://sbml.org/Main_Page) models.
Performs the following checks:

- Use of parameters or species that are uninitialized (and
not present on the right hand side of a reaction)
- Correctly constructed kinetics laws
- Parameters that are unreferenced
- Mass conservation
