PocketOptimizer - A Python Library for Protein-Ligand Binding Design
====================================================================

PocketOptimizer is a framework that offers experimentation with different scoring functions in
protein-ligand binding design.

- **API** -- Python 3 and workflow interface
- **Minimization** -- GPU based minimization through `OpenMM <https://openmm.org/>`_
- **Force Fields** -- `Amber ff14SB <https://pubs.acs.org/doi/10.1021/acs.jctc.5b00255>`_ and `CHARMM36 <https://pubmed.ncbi.nlm.nih.gov/23832629/>`_
- **Scoring Functions** -- Binding interaction scoring with `Smina <https://github.com/mwojcikowski/smina>`_ (empirical) or `FFEvaluate <https://software.acellera.com/htmd/tutorials/FFEvaluate.html>`_ (physics-based)
- **Deterministic** -- Linear programming is applied to find the Global Minimum Energy Conformation (GMEC) in a given set

Installation via the Repository
-------------------------------

If you want to use PocketOptimizer you can clone the GitHub repository and set up the conda environment:

.. code-block:: bash

    git clone https://github.com/Hoecker-Lab/pocketoptimizer

    cd /YOUR_PATH/pocketoptimizer

The repository contains an environment.yml file that holds information about all
dependencies and versions required for PocketOptimizer. By typing:

.. code-block:: bash

  conda env create -f environment.yml

Conda creates a separate PocketOptimizer environment.

Alternatively, you can use the minimal environment:

.. code-block:: bash

  conda env create -f simple_env.yml

**Warning**: If you want to use CUDA, you might encounter a compatibility problem because the
CUDA version in the environment.yml is not compatible with your driver. Therefore use the
simple_env.yml instead.

Running PocketOptimizer
=======================

Since PocketOptimizer's rework, the framework is now accessed solely through Python 3.9.
This means that it can be incorporated in your regular Python scripts or in an interactive fashion through Jupyter Notebooks.
The former is useful for automatized pipelines while the latter is helpful for going through every step separately.
Jupyter is already part of the ``pocketoptimizer`` environment and does not need be installed.
Jupyter can be opened in two different ways:

``jupyter-lab``

or

``jupyter-notebook``


Before being able to use PocketOptimizer in Jupyter you need to install the kernel from inside the conda environment:

.. code-block:: bash

    conda activate pocketoptimizer

    python3 -m ipykernel install --user --name pocketoptimizer

Then you can select it in the upper right corner of the Jupyter-Notebook.

General Project Layout
----------------------

PocketOptimizer works for each protein-ligand design project in a separate project
directory.

The general layout looks like this:

::

    project
    ├── designs
    ├── energies
    ├── ligand
    └── scaffold

The first step involves creating a ``project`` directory (you can obviously choose the name)
and a ``scaffold`` and ``ligand`` sub-directory.
The other directories are automatically created during the design process.

Now place the following files inside the ``ligand`` and ``scaffold`` directory:

::

    project
    ├── ligand
    │   └──  YOUR_LIGAND.mol2
    │
    └── scaffold
        └── YOUR_PROTEIN.pdb

* ``YOUR_LIGAND.mol2`` = starting ligand pose placed inside the binding pocket
* ``YOUR_PROTEIN.pdb`` = protein structure used as scaffold

1. Ligand Preparation
---------------------

1.1 How to get your small molecule
**********************************

There are multiple ways to obtain your molecule of choice.
If you want to make a design for a molecule different from
a ligand bound in your crystal structure, you can do a search on
`RCSB <http://www.rcsb.org/pdb/ligand/chemAdvSearch.do>`_ for different kinds of ligands.
This allows you to download a molecule in the sdf format.

If you already have a protein crystal structure with the desired ligand, you can also
extract the ligand from the .pdb file using for example `PyMol <https://pymol.org/2/>`_. But beware that the ligand
is missing all hydrogen atoms.

**Note**: PocketOptimizer works with several input formats (mol2, sdf) that will be converted internally.


1.2 Placing the ligand inside the binding pocket
************************************************

PocketOptimizer is based on semi-rational design principles which offers the
flexibility to design the binding pocket following your ideas.

If you extracted your ligand from a protein crystal structure, then this step is
not of importance for you. Otherwise, the easiest way to get the ligand inside the binding pocket is to superpose it
on an existing ligand. The superposition is strictly dependent on your design
thoughts and also requires some experimentation und multiple design runs.

The easiest way the superposition can be done is to use PyMol, which offers
a Pair-Wise alignment tool to easily align elements the way you want to. The tool
can be found in the PyMol toolbar at the top in ``Wizard`` as the
name ``Pair Fit``.

If you don't have initial information about a binding pose available, another way is to produce an initial
pose using a docking program such as `Autodock Vina
<https://vina.scripps.edu/>`_.


2. First Design Steps
---------------------

As mentioned, PocketOptimizer needs to be initialized in your project directory.
Therefore, inside every script or Jupyter notebook you use, you need to define
the following lines:

.. code-block:: python

    # Append the PocketOptimizer Code
    import sys
    sys.path.append('YOUR_POCKETOPTIMIZER_PATH')

    # Import the pocketoptimizer module
    import pocketoptimizer as po

    # Initialize a new design pipeline
    design = po.DesignPipeline(work_dir=project_dir,         # Path to working directory containing scaffold and ligand subdirectory
                               ph=7,                         # pH used for protein and ligand protonation
                               forcefield='amber_ff14SB',    # forcefield used for all energy computations (Use Amber as it is better tested!)
                               ncpus=8)                      # Number of CPUs for multiprocessing

While you are initializing you can define a ``pH``, used for protonating the side chains of the protein and also the ligand molecule.
Additionally, PocketOptimizer has two ``force fields`` implemented, the AMBER ff14SB and the CHARMM36 force field.
These force fields contain parameters and energy functions to calculate the energy of the protein-ligand system.
Besides you can define the ``number of CPUs`` used for all energy calculations.

2.1 Preparation/Minimization
****************************

2.1.1 Ligand Preparation
++++++++++++++++++++++++

The ligand also gets protonated and parameterized. However, the chemical space for small molecules
can not be easily described by prebuild force field atom types, since the variety of small organic
molecules far exceeds that of the 20 canonical amino acids, which is why ligands generally need to be
parameterized separately. For AMBER force fields this can be done by using either `GAFF or GAFF2 (General
AMBER Force Field) <https://pubmed.ncbi.nlm.nih.gov/15116359/>`_, for CHARMM the tool is called
`CGenFF (Charmm GENeral Force Field) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2888302/>`_.

PocketOptimizer needs the following ligand inputs:

    * Ligand in `mol2 <https://zhanggroup.org//DockRMSD/mol2.pdf>`_/`sdf <https://chem.libretexts.org/Courses/Intercollegiate_Courses/Cheminformatics_OLCC_(2019)/2._Representing_Small_Molecules_on_Computers/2.5%3A_Structural_Data_Files>`_ format

Eventually:

    * Parameters in `frcmod <https://ambermd.org/FileFormats.php#frcmod>`_ or `prm <https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node25.html>`_/`rtf <https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node24.html>`_ format

Experienced users can obtain these by using tools like `ANTECHAMBER <http://ambermd.org/antechamber/ac.html>`_ and
`PARMCHK <http://ambermd.org/tutorials/basic/tutorial5/>`_ for the AMBER force field or `CGenFF <https://cgenff.umaryland.edu/>`_ for the CHARMM force field.

PocketOptimizer offers a python interface utilizing these tools to parameterize your small molecule:

.. code-block:: python

    #  Only necessary if you don't have ligand parameters.
    design.parameterize_ligand(
    input_ligand='ligand/YOUR_LIGAND.mol2', # Input ligand structure file could be .mol2/.sdf
    addHs=True                              # Whether to add hydrogen atoms to the input structure
    )

This creates a ``ligand.mol2`` structure file and additionally either a ``ligand.frcmod`` or ``ligand.prm``/``ligand.rtf`` parameter files in the ``ligand``
directory under ``FORCE_FIELD/params``. Before you proceed, take a look at those files if the structure is correct protonated and suits your needs.

::

   ligand
   ├── ligand_structure.mol2
   └── FORCE_FIELD
       ├── ligand.mol2
       └── params
           └── ligand.mol2/ligand.frcmod or ligand.prm/ligand.rtf

**Hint**: Use relative paths for the scaffold and ligand structures,
as you are inside the project directory during the entire design process.

2.1.2 Protein Preparation
+++++++++++++++++++++++++

Before the design process can start, the protein scaffold needs to be cleaned of ions, waters, small molecules (like natural ligands)
and unnecessary protein chains. Furthermore, the protein scaffold needs to be protonated to a certain pH that was defined when
initializing the design pipeline and it needs to be minimised. This is because experimentally solved protein structures commonly
do not contain hydrogen atoms and often have clashes due to crystallographic model building.
PocketOpimizer has built in functionalities for this, utilizing the `HTMD <https://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00049>`_
and `OpenMM <https://openmm.org/>`_ distribution. After you placed your protein of choice inside the ``PROJECT_NAME/scaffold/``
directory, you can start to open a Python/IPython console or preferably a Jupyter
notebook and type the following:

.. code-block:: python

    design.prepare_protein(
        protein_structure='scaffold/YOUR_PROTEIN.pdb',  # Input PDB
        keep_chains=['A', 'B'],  # Specific protein chain to keep
        backbone_restraint=True, # Restrains the backbone during the minimization
        cuda=False,              # Performs minimization on CPU instead of GPU
        discard_mols=[]          # Special molecules to exclude. Per default everything, but peptides have to be defined manually
        )

This allows to minimize the structure with or without the backbone being constrained.
Remember, this can also be a design choice you want to consider
as the scaffold/backbone is the foundation of your design.
The following files are created after this step:

::

    scaffold
    └── FORCE_FIELD
        ├── protein_preparation
        │   ├── prepared_scaffold.pdb
        │   └── scaffold_report.xlsx
        ├── protein_params
        └── scaffold.pdb

In the scaffold folder a ``FORCE_FIELD`` sub-folder is created named after the respective
force field that was set in the beginning of the design process. Within this folder, a
``protein_preparation`` sub-folder is created, which contains the cleaned and protonated protein structure.
A scaffold report in form of an excel spreadsheet is also created within this folder that
contains information about the modified residues (like protonation states or filled-in missing atoms (hydrogen atoms)).

A ``protein_params`` sub-folder is created within the ``FORCE_FIELD`` sub-folder that contains force field parameters and energy
functions describing the protein, which can be used to calculate various interaction-energies.

2.2 Choose your design positions
********************************

Next you can start taking a look at the resulting structure in:

::

    scaffold
    └── FORCE_FIELD
        └── scaffold.pdb


This is the protonated and minimized version of your initial protein, you can start to choose the
residues you want to mutate or you want to be flexible:

.. code-block:: python

    # Your mutations
    design.set_mutations([
        {'mutations': ['ALA', 'ASN', 'GLU'], 'resid': '8', 'chain': 'A'},
        {'mutations': ['LEU'], 'resid': '10', 'chain': 'A'},
        {'mutations': ['SER'], 'resid': '12', 'chain': 'A'},
        {'mutations': ['TYR'], 'resid': '28', 'chain': 'A'},
        {'mutations': ['PHE'], 'resid': '115', 'chain': 'A'},
    ])

The design positions are defined as a list containing dictionaries for every
design position. If only a single amino acid is provided in the mutations list, only a single
option is tested. This can be used to model the flexibility of native residues
you don't want to mutate, but instead to move (rotate). Residues not defined
in this list are static during the design and don't move at all.

You can also use certain keywords to try out a number of amino acids, grouped by their properties:

.. code-block:: python

        'ALL': ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HID', 'HIE', 'HIP',
                'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'],
        'AROMATIC': ['PHE', 'TRP', 'TYR'],
        'AMIDE': ['ASN', 'GLN'],
        'ALIPHATIC': ['GLY', 'ALA', 'VAL', 'LEU', 'ILE'],
        'ACIDIC': ['ASP', 'GLU'],
        'BASIC': ['LYS', 'ARG'],
        'HYDRO': ['SER', 'THR'],
        'SULF': ['CYS', 'MET']

Once you are done and the mutations are defined, you can start preparing the
mutated scaffolds for the later energy and scoring calculations
(the parameters for the prepared scaffolds are also contained within
the ``protein_params`` sub-folder):

.. code-block:: python

    # Prepares all defined mutants and glycine scaffolds for side chain rotamer and ligand pose sampling
    design.prepare_mutants(sampling_pocket='GLY')

**Hint**: Testing additional residues/mutations later on is not a problem.
PocketOptimizer dynamically detects which mutations/calculations already exist and only calculates additional ones.

**Note**: If you add or remove design positions, you will need to create an entirely new design.

3. Sampling Flexibility
-----------------------

The following steps are definitely the most time consuming ones and have therefore
the option to be multiprocessed.

The steps that are now needed contain:

* Calculation of possible rotamers
* Calculation of possible ligand poses
* Computation of the energies and scores

3.1 Create Ligand Conformers
****************************

To model your ligands flexibility correctly, a .pdb file containing ligand conformations is needed.

::

     ligand
     └── FORCE_FIELD
         └── conformers
             └── ligand_confs.pdb

PocketOptimizer has an interface for `Obabels <https://open-babel.readthedocs.io/en/latest/3DStructureGen/multipleconformers.html>`_ conformer sampling:

.. code-block:: python

        # Obabel conformer generation
        design.prepare_lig_conformers(
        nconfs=50,         # Maximum number of conformers to produce (Sometimes these methods produce lower number of conformations)
        method='genetic',  # Genetic method in OpenBabel, other option is confab
        score='rmsd',      # For genetic method: filters conformers based on RMSD diversity or filtering based on energy diversity
        #rcutoff=0.5,  # Confab method: RMSD cutoff
        #ecutoff=50.0 # Confab method: Energy cutoff
        )


This samples a maximum number of 50 conformers using either a ``genetic`` algorithm or
the ``confab`` procedure as implemented in Obabel. The ``genetic`` algorithm derives
at an optimal solution either based on RMSD or energy diversity after a series of generations.
The ``confab`` method systematically generates conformers based on a set of allowed torsion angles
for every rotatable bond and prunes out conformers based on an energy threshold and RMSD diversity.

3.2 Create Ligand Poses
***********************

The ligand pose sampling procedure requires the user to define a grid that specifies
in which range possible ligand poses are going to be sampled. This procedure generates a number of poses from the
sampled ligand conformers by translating and rotating them along a user defined grid. Alternatively, a random sampling procedure can
be performed by setting the parameter method to ``random``.

.. code-block:: python

    # Sampling of ligand poses
    # Defines a grid in which the ligand is translated and rotated along.
    #                       Range, Steps
    sample_grid = {'trans': [1, 0.5],  # Angstrom
                   'rot': [20, 20]}    # Degree
    design.sample_lig_poses(
        method='grid',         #  Uses the 'grid' method. Other option is 'random'
        grid=sample_grid,      #  Defined grid for sampling
        vdw_filter_thresh=100, #  Energy threshold of 100 kcal/mol
        max_poses=10000        #  Maximum number of poses
    )

The grid is defined in a Python dictionary that containes rotational and translational
movements in the following form ``[MAXIMUM DISTANCE/ANGLE, STEPS]``, which means
that in the shown example the ligand would be moved 1 angstrom around every axis
in 0.5 angstrom steps and rotated by 20 degree around every axis in 20 degree steps.
A vdW energy threshold ensures that the sampled poses are not clashing with the
scaffold. This ligand pose pruning procedure is again performed in a glycine scaffold,
where all design positions are mutated to the amino acid glycine. If the number of
accepted poses exceeds the maximum number of poses defined, a MinMax diversity Picker
from RDKit will be applied to filter all sampled poses based on maximum RMSD diversity.

The ligand poses are saved as frames of a trajectory in the files ``ligand_poses.pdb``
and ``ligand_poses.xtc``. Furthermore, their energies can be inspected in ``ligand_poses.csv`` under:

::

     ligand
     └── FORCE_FIELD
         └── poses
             ├── ligand_poses.pdb
             ├── ligand_poses.xtc
             └── ligand_poses.csv


3.3 Create Side Chain Conformers
********************************

Side chain rotamers can be sampled with the following method based on the fixed backbone that has been prepared:

.. code-block:: python

    # Sampling of side chain rotamers
    design.sample_sidechain_rotamers(
        vdw_filter_thresh=100,         # Energy threshold of 100 kcal/mol for filtering rotamers
        library='dunbrack',            # Use dunbrack rotamer library (Should be used!)
        dunbrack_filter_thresh=0.01    # Probability threshold for filtering rotamers (1%)
        )

This procedures will use the design mutations that were set in the previous step and a defined van
der Waals energy threshold to prune rotamers that clash with the protein scaffold.
The default value is 100 kcal/mol. This pruning procedures are
also performed in your defined sampling scaffold (glycine), where all other design positions are
mutated to the amino acid glycine.

Additionally, a rotamer library can be selected.
Options are either the original PocketOptimizer rotamer library ``CMLib`` or the backbone dependent
`Dunbrack rotamer library <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3118414/>`_.
When using the Dunbrack rotamer library a filter threshold can be defined which allows
to filter out all rotamers that have a probability of occuring of less than the defined threshold.
Accordingly, the threshold should be between 0 and 1 and allows to reduce the amount of sampled rotamers.

All accepted rotamers are contained in .pdb files and their energies are contained in .csv files under:

::

    scaffold
    └── FORCE_FIELD
        ├── scaffold.pdb
        └── rotamers
            └──  LIBRARY
                 └── POSITION
                     ├── RESNAME.csv
                     └── RESNAME.pdb

4. Energy Calculations
----------------------

Next all protein-protein and protein-ligand interaction energies are calculated, the protein-protein interaction energies are evaluated from force fields,
whereas the protein-ligand interaction energies can be also evaluated using different scoring functions. To get an overview over all available scoring functions:

.. code-block:: python

    # Outputs all available scoring functions
    design.scoring
    {'smina': ['vina', 'vinardo', 'ad4_scoring'],
     'ff': ['amber_ff14SB', 'charmm36']}

To calculate the energies:

.. code-block:: python

    # Calculate the binding and packing energies of all ligand poses and side chain rotamers against each other and against the fixed scaffold
    design.calculate_energies(
        scoring='vina', #  Method to score protein-ligand interaction
        )

This step also defines the used scoring function (to change it repeat the step and use a different scoring function).
All energies are contained in .csv files under:

::

    project
    ├── designs
    ├── energies
    │   └── FORCEFIELD_LIBRARY
    │       ├── sidechain_scaffold_FORCE_FIELD
    │       │   └── RESIDUE.csv
    │       ├── sidechain_sidechain_FORCE_FIELD
    │       │   └── RESIDUE_A_RESIDUE_B.csv
    │       ├── ligand_scaffold_SCORING_METHOD
    │       │   └── ligand.csv
    │       └── ligand_sidechain_SCORING_METHOD
    │           └── ligand_RESIDUE_A.csv
    ├── ligand
    └── scaffold


5. Design Solutions
-------------------

After the energy computations are finished, the best ligand poses/rotamers can be
calculated in order to finish the PocketOptimizer run.

This is where PocketOptimizer shines the most, because you have a lot of freedom
to experiment with the force field and scoring functions you used before and also
how to scale them.

The final designs can be calculated with:


.. code-block:: python

    # Compute the lowest energy structures using linear programming
    design.design(
        num_solutions=10,           #  Number of solutions to compute
        ligand_scaling=10,          #  Scaling factor for binding-related energies (You need to adapt this to approximate the packing and binding energies)
    )

which first prepares input files for the optimizer and then creates output
.html/.txt files and pymol sessions containing all the designed structures:

::

    project
    ├── designs
    │   └── FORCE_FIELD_SAMPLING_LIBRARY
    │       └── DESIGN_MUTATIONS
    │           └── SCORING_METHOD_LIGAND_SCALING
    │               ├──  INDEX_DESIGN_SOLUTION
    │               │    ├── ligand.mol2
    │               │    ├── receptor.pdb
    │               │    ├── report.txt
    │               │    ├── report.html
    │               │    └── design.pml
    │               ├── summary.txt
    │               ├── summary.html
    │               ├── summary.pml
    │               ├── summary.png
    │               └── seqlogo.png
    ├── energies
    ├── ligand
    └── scaffold

Every design solution is contained as a single folder named after the index of the solution,
this folder contains a structure for the receptor and ligand of the design respectively as
well as the reports and a pymol session. Summaries of the energies for all best design solutions
are contained in summary.txt/.html files and all the structures are contained in a summary pymol
session. All energies are also graphically depicted in a summary energy plot. If multiple residues are
allowed at design positions, a sequence logo is generated. The sequence logo depicts
design position together with the frequency of mutations at these positions.

**Note**: It is important to take a look at the energies contained in the .txt/.html and
also to inspect the final output structures.

5.2 Cleaning the working directory
**********************************

PocketOptimizer creates many files in the directory that is specified as the working directory.
These can be files containing parameters for the protein or the ligand molecule or files containing the calculated energies.
In order to delete them, PocketOptimizer includes a clean-up procedure, which scans your working directory after these files.

.. code-block:: python

    design.clean(
        scaffold=True, #  Deletes all scaffold-related files
        ligand=True    #  Deletes all ligand-related files
    )

You can specify if you want to delete only the files related to the scaffold or the ligand or both. This deletes all files
that were created during the design run and allows you to start an entirely new design in your working directory.


6. Command Line Interface
-------------------------

By running the python script: ui.py, you can also access the command line interface:

.. code-block:: bash

    usage: ui.py [-h] [-ff FORCEFIELD] -r RECEPTOR -l LIGAND [--ph PH] --mutations MUTATIONS [MUTATIONS ...] [--vdw_thresh VDW_THRESH] [--library LIBRARY]
                 [--nconfs NCONFS] [--rot ROT] [--rot_steps ROT_STEPS] [--trans TRANS] [--trans_steps TRANS_STEPS] [--max_poses MAX_POSES]
                 [--scoring SCORING] [--scaling SCALING] [--num_solutions NUM_SOLUTIONS] [--ncpus NCPUS] [--cuda] [--clean]

    PocketOptimizer CLI, for more options use API.

    optional arguments:
      -h, --help            show this help message and exit
      -ff FORCEFIELD, --forcefield FORCEFIELD
                            Force field to be used either: amber_ff14SB or charmm36
      -r RECEPTOR, --receptor RECEPTOR
                            Protein input structure file in pdb format
      -l LIGAND, --ligand LIGAND
                            Ligand input structure file
      --ph PH               ph value for side chain and ligand protonation
      --mutations MUTATIONS [MUTATIONS ...]
                            Mutations (A:1:ALA)
      --vdw_thresh VDW_THRESH
                            Energy threshold for rotamer and ligand pose sampling (kcal/mol)
      --library LIBRARY     Rotamer library, options are: dunbrack or cmlib
      --nconfs NCONFS       Number of ligand conformers to sample, default: 50
      --rot ROT, --rot ROT  Maximum ligand rotation, default: 20°
      --rot_steps ROT_STEPS, --rot_steps ROT_STEPS
                            Ligand rotation steps, default: 20°
      --trans TRANS, --trans TRANS
                            Maximum ligand translation, default: 1 Å
      --trans_steps TRANS_STEPS, --trans_steps TRANS_STEPS
                            Ligand translation steps, default 0.5 Å
      --max_poses MAX_POSES, --max_poses MAX_POSES
                            Maximum number of ligand poses to sample, default: 10000
      --scoring SCORING     Scoring function, options are: vina, vinardo, ad4_scoring, amber_ff14SB or charmm36
      --scaling SCALING     Ligand scaling factor, default: 1
      --num_solutions NUM_SOLUTIONS
                            Number of design solutions to calculate, default 10
      --ncpus NCPUS         Number of CPUs for multiprocessing
      --cuda                Enabling cuda for GPU-based minimization
      --clean               Clean the working directory



LICENSE
-------

PocketOptimizer is licensed under the GNU GENERAL PUBLIC LICENSE however we would like to point out
the `HTMD Software Academic License Agreement <https://github.com/Acellera/moleculekit/blob/main/LICENSE>`_
which makes the software for academic use only.


Publications
************

**PocketOptimizer 2.0: A modular framework for computer-aided ligand-binding design**, Noske J, Kynast JP, Lemm D, Schmidt S, Höcker B.
Protein Sci. 2022 Nov 19:e4516. doi: `10.1002/pro.4516 <https://pubmed.ncbi.nlm.nih.gov/36403089/>`_.

**PocketOptimizer and the Design of Ligand Binding Sites**, Stiel AC, Nellen M, Höcker B.,
Methods Mol Biol. 2016;1414:63-75. doi: `10.1007/978-1-4939-3569-7_5
<https://www.ncbi.nlm.nih.gov/pubmed/27094286>`_.

**Binding pocket optimization by computational protein design**, Malisi C, Schumann M, Toussaint NC, Kageyama J, Kohlbacher O, Höcker B.,
PLoS One. 2012;7(12):e52505. doi: `10.1371/journal.pone.0052505
<https://www.ncbi.nlm.nih.gov/pubmed/23300688>`_.
