# disp_solve

disp_solve is a script that generates a number of atomic displacement cell files for a crystal compound. These
displaced unit cell atoms result from enacting a range of uniform electric fields on a system. It takes input files
from a CRYSTAL simulation to generate the necessary matrices for this displacement calculation.

## Installing python 3

A **python 3** distribution must be installed. Either Anaconda or direct installation [here](https://realpython.com/installing-python/) are the best options. 

## Correct packages

The following packages must be installed with your distribution for the program to run:

- numpy
- sys
- itertools
- copy
- os.path
- tabulate
- ast
- pathlib
- json

Install packages you don't have with `pip install <module_name>` or check that you have it with `pip list | grep <module_name_you_want_to_check>`

## Input files

The program requires the following inputs:

- CRYSTAL output file of a phonon calculation of your system including IR intensities(Born charges) and Hessian matrix printed in fractional coordinates
- The initial system.cell file of the ENV code corresponding to the undistorted system 

The program has optional input:

- disp_solve.input command line input file containing the necessary data from the user (e.g file locations, electric field ranges, input method and charge data if needed)

## Executing script

There are two options:

1. `python3 disp_solve.py` (no input file: command line will prompt you for direct inputs. Can use full or relative path for script)
2. `python3 disp solve.py <input_file_path>` (Input file is passed as a command line argument. Can be relative or full path.)

If you choose 1. you will be prompted for 3-4 different things:

1. Unit cell generator parameter: direct, charge or auto
   - Entering `direct` makes the program extract the unit cell irreducible group atom information, unit cell coordinates and charge data from the same initial system.cell file.
   - Entering `auto` generates all of the unit cell information from the CRYSTAL output file, except the charges which come from the system.cell file (useful if this initial cell file is in a different spatial group).
   - Entering `charge` does the same as `auto` but prompts manual user input of the charge for each individual atom in the irreducible group.
2. Relative/full path directories of the CRYSTAL output file and the initial system.cell file
3. Start,stop,step of the ranges for E\_a, E\_b, E\_c inputs to generate the uniform electric field cell grid (given in kV/m).
4. If the `charge` option is chosen, charges are prompted for. 

If you choose 2. you need to write the correct format for the third input file. Here is an example:

```
crystal_file = ht.frequence.B1PW_PtBs.loto.out
cell_init = ymno3.cell
ea = [0,0.5,0.1]
eb = [0,0.2,0.05]
ec = [0,1,0.2]
unit_source = charge
charge_dict = {"Y":3.0, "MN":3.0, "O1":-2.0, "O2":-2.0}
```

- Each line element is separated by a space (the = sign is isolated)
- Crystal and cell file relative or full path directories are assigned to crystal\_file and cell\_init respectively
- ea, eb, ec are the ranges for E\_a, E\_b, E\_c respectively in the format \[start,stop,step\]
- unit_source is assigned to either direct, auto or charge
- If charge is chosen above, give each atom's charge in the [python dictionary format](https://www.w3schools.com/python/python_dictionaries.asp)

## Output

- The output file is named after the initial cell file with an added (Ea,Eb,Ec) component label.
- If the files are call locally, output cell file grid will be stored in the local directory folder named after the array parameters for E\_a, E\_b, E\_c.
- If the full path is called, output grid will be generated in the same way but in the provided full path directory
