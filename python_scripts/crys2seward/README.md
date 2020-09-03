# Crys2Seward

Crys2seward is a python program that links the output from the CRYSTAL geometry optimisation to `disp_solve`, `env15` and `env2seward` to produce a grid of input files for the program SEWARD.

## Installing python 3

A **python 3** distribution must be installed. Either Anaconda or direct installation [here](https://realpython.com/installing-python/) are arguably the best options. 

## Correct packages

The following packages must be installed with your distribution for the program to run:

- sys
- os
- shutil
- pathlib
- env2seward
- disp\_solve

Install packages you don't have with `pip install <module_name>` or check that you have it with `pip list | grep <module_name_you_want_to_check>`, unless they are packages available only in this repository such as env2seward and disp\_solve.

In that case you will need to download both the script and the `__init__.py` files, and reference their path in crys2seward.py when declaring used packages (unless they are located in the same relative path as the repository, then no changes are required).

## Input files

The program requires the following inputs:

- CRYSTAL output file of a phonon calculation of your system including IR intensities (Born charges) and Hessian matrix printed in fractional coordinates
- The initial system.cell file of the ENV code corresponding to the undistorted system
- The envin file for the env15 program `xenv15` after executing the necessary `make` command in the Env15 folder (see the local tools manual written by M. B. Lepetit for how to install env15).
- The prefix.c2s.in input file containing all the paths pointing to the aforementioned files, as well as the inputs for `disp_solve` and `env2seward`

## Executing the script

To execute the script, type:

`python3 <path_to_crys2seward.py> <path_to_c2s.in_input_file>`

into the command line in the same working directory as the xenv15 program and it's input files (envin and env.out).
xenv15 only works in the cwd, i.e. it cannot be referenced from another path.

For example, on this git repository the relative paths are:

`python3 ../../../python_scripts/crys2seward/crys2seward.py ../ymno3_d1.c2s.in`

The main difficulty of executing the script is writing a correct input file. Here is an example:

### ymno3\_d1.c2s.in

```
#########crystal2seward input########################
envin = ymno3_d1.envin
envout = ymno3_d1.env.out
xenv_sew0 = ymno3_d1.sew0
xenv_psd = ymno3_d1.psd
#########disp_solve input############################
# Born and Hessian sources can either be both '.DAT' or both '.loto.out' filetype
born_file = ../ex_student_ymno3_data/Position/YMnO3/BORN_B1Pw_loto.DAT
hess_file = ../ex_student_ymno3_data/Position/YMnO3/HESSIEN.DAT
crystal_file = ../ex_student_ymno3_data/frequence.B1PW.loto.out
cell_init = ../ymno3.cell
ea = [0,.15,.03]
eb = [0,.15,.03]
ec = [0,.15,.03]
unit_source = auto
charge_dict = {'Y':3.0, 'MN':3.0, 'O1':-2.0, 'O2':-2.0}
##########env2seward input###########################
filename = ymno3_d1.0.0_0.0_0.0.sew.in
title = <insert_title_here>
sew0_file = ymno3_d1.0.0_0.0_0.0.env.sew0
psd_file = ymno3_d1.0.0_0.0_0.0.env.psd
lib_frag = {'O': {'key': 'O_frag_basis', 'loc': ''}, 'Y': {'key': 'Y_frag_basis', 'loc': ''}, 'MN': {'key': 'MN_frag_basis', 'loc': ''}}
lib_pseudo = {'MN': {'key': 'MN_pseudo_basis', 'loc': 'MN_pseudo_library'}, 'O1': {'key': 'O1_pseudo_basis', 'loc': 'O1_pseudo_library'}, 'O2': {'key': 'O2_pseudo_basis', 'loc': 'O2_pseudo_library'}, 'O3': {'key': 'O3_pseudo_basis', 'loc': 'O3_pseudo_library'}, 'O4': {'key': 'O4_pseudo_basis', 'loc': 'O4_pseudo_library'}, 'Y': {'key': 'Y_pseudo_basis', 'loc': 'Y_pseudo_library'}, 'Y1': {'key': 'Y1_pseudo_basis', 'loc': 'Y1_pseudo_library'}, 'Y2': {'key': 'Y2_pseudo_basis', 'loc': 'Y2_pseudo_library'}}
```

### Rules for c2s.in:

- envin and envout correspond to the absolute path of the prefix.envin, prefix.env.out files
- xenv\_sew0 and xenv\_psd correspond to the default names of the xenv15 output files for the system
- The disp\_solve and env2seward input rules can be found in their README's and the software manual

## Output

Each grid of files is stored in an appropriately named folder. The sew0, psd and sew.in files are generated in separate folders but in the same current working directory with the xenv15 files as discussed before.
The disp\_solve cell files are generated in the same way, but stored in the same directory as the initial cell file or the CRYSTAL output file if their relative or full paths are specified.
They are **not** generated in the local directory unless the disp\_solve input file paths are not specified. Every grid folder name features the start, stop, step for Ea, Eb, Ec chosen in the c2s.in file in the disp\_solve input section.

## How does it work?

The CRYSTAL optimised geometry output files contains the associated minimum energy Hessian and Born charge tensor data, alongside the unit cell atomic coordinates and lattice parameters.
CRYSTAL also produces the unit cell coordinates alongside the atom's charges in a **cell file**.

This data is used to calculate the displacements from an induced uniform electric field by `disp_solve`, and a grid of new cell files are generated with electrically disturbed positions, based on the initial cell file passed as an input.

Crys2seward edits the xenv15 'envin' or 'env.in' input files so that the env15 program can take cell files of different displacements and produces a grid of output files named `prefix.env.sew0` and `prefix.env.psd`.
These are the input files of the env2seward program, where the prefix corresponds to the simulated system.

Crys2seward modifies the env2seward input data so that each set of sew0 and psd file coordinates are categorised into a prefix.sew.in format.

Each file grid (cell files, env15 output files, sew.in files) is contained in a directory named after the electric field ranges and filetype. The names of each grid file include the associated incident electric field.


