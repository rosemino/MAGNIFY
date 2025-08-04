This directory contains the force field, well-tempered metadynamics parameters, LAMMPS input files and equilibrated configurations for the [Zn(MIm)₂] system used in the well-tempered metadynamics runs, as well as the data file for the unbiased nucleation systems filled with MeOH or DMSO (reference system) from the publication: "Unveiling ZIF-8 nucleation mechanisms through molecular simulation: role of temperature, solvent and reactant concentration" by Sahar Andarzi Gargari and Rocio Semino

Atom types:
 type 1 # Zn
 type 2 # C2
 type 3 # N1
 type 4 # H2
 type 5 # H3
 type 6 # C1
 type 7 # C3
 type 8 # He
 type 9 # Ne

He and Ne are dummy atoms attached to the Zn and N species to correctly reproduce the angular distribution of ligands around a Zn center.

Solvent atom types:
DMSO:
 type 10 # S4
 type 11 # O4
 type 12 # C4 (CH₃ group defined under a single name: C4)

MeOH:
 type 10 # C (CH₃ group defined under a single name: C)
 type 11 # O
 type 12 # H

Walkers:
 To run multiple walkers, separate directories must be prepared, and the WALKERS_ID must be set uniquely in each PLUMED input file; for example, with 5 walkers, WALKERS_ID should be defined as 0, 1, 2, 3 and 4 across the directories.

Reference:
http://arxiv.org/abs/2507.23574
