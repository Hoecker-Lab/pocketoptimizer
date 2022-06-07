# Preparation Pipeline for PocketOptimizer

This tools here are not necessary if you already have a minimized protein structure
built with hydrogens, protonations and ligand conformers in .pdb format.

Otherwise this tools here can guide you on the way. However, further calculations in
PocketOptimizer require certain structure parameters that are generated in the structure
building procedure.

## structure_building

Contains tools for side chain protonation and structure building with Amber/Charmm.
Does also contain a short ligand parameterization interface.

## minimize_structure

Structure minimization.

## conformer_generator_obabel

Conformer generation with OpenBabel using either the genetic or confab method.

## aacodes

Mapping between non-canonical and canonical amino acids.