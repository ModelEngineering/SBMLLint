# SbChecker
Model checker for the Antimony kinetics modeling language.
Performs the following checks:

- Use of parameters or species that are unitialized (and
not present on the right hand side of a reaction)
- Correctly constructed kinetics laws
- Parameters that are unreferenced
- Mass conservation
