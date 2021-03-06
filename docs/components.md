# Component Specifications

<img src="components.png" width="200"/>

This document describes the major classes in SBMLLink and the relationships between them.

The diagram depicts the component hierarchy.
If component A is above B then A knows about B but but does not know about A.
Components at the same level do not know about each other.

Below, we describe each component.

- Moiety. A chemical group within a molecule.

- MoietyStoichiometry. A Moiety with a count of its occurrence.

- Molecule. A representation for a chemical structure. A molecule has a name that uniquely identifies the object.

- MoleculeStoichiometry. A Molecule with its count of occurrence within a collection.

- Reaction. A representation of the transformation of a set of molecules (reactants) into another set of molecules (products).
A reaction is uniquely identified by an identifier (a string representation
of the reaction, including its kinetics law).

- SimpleModel. A representation of the SBML model.

- SOM. A set of molecules.

- MoietyComparator. Reports on the differences between two occurrences of moieties within a collection of MoleculeStoichiometry.

- Arc. Indicates a mass inequality relationship between two SOMs.

- MESGraph. Used to make inferences about the inequalities of masses of molecules.

- Tool. A command line analysis done on an SBML file.

We refine the foregoing by looking at the more detailed relationships
between classes in the following diagram.

<img src="reaction_molecule_moiety_class_diagram.png" width="200"/>

The graph can be summarized as follows.

A Reaction has a kinetics, law, a unique identifier, Reactants,
and Products. The latter two are collections of MoleculeStoichiometry,
a molecule with its stoichiometry.
A Molecule has a name. It can be structured
into one or more MoietyStoichiometry. 
For example, the molecule `ATP` can be written
in terms of its moieties as `A_P_P_P` or, using repetition counts,
`A__P_3`. `A` and `P` are moieties, and `P_3` is a moiety with its
stoichiometry.
