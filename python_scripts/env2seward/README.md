# Env2seward
This code processes and formats the output from ENV so it can be used in a calculation
performed by SEWARD.

## Installing python 3

A **python 3** distribution must be installed. Either Anaconda or direct installation [here](https://realpython.com/installing-python/) are the best options. 

## Input files

The program requires the following inputs:

- prefix.env.sew0 is a help file that contains the input data for SEWARD of the MOLCAS chain. It contains the data for the quantum fragment, the TIPS and the renormalised charges
- prefix.env.psd contains the information on the atoms represented by TIPs

The program has optional input:

- prefix.input command line input file containing the necessary data from the user (e.g file paths, titles, basis set libraries)

## Executing the script

There are two options:

1. `python3 env2seward.py` (no input file: command line will prompt you for direct inputs. Can use full or relative path for script.)
2. `python3 env2seward.py <input_file_path>` (Input file is passed as a command line argument. Can be relative or full path.)

For option 1. you will be prompted for:

1. Name of the system prefix such as `GdMn2O5_J1' (Coordinate input files must be local for this option)
2. Title and name of created file
3. Basis sets and library locations

If you want to call the input files from full path, or automate the process in a batch/ shell file, use option 2.

For option 2 you will need to write the third input file with the following format:
```
filename = chosen_filename
title = chosen_title
# Following files can be local path or full path, just replace name with full path directory
sew0_file = prefix.env.sew0
psd_file = prefix.env.psd
lib_frag = {"atom type":{"loc":"specified basis set library","key":"basis set for
atom type"},"atom type2":...}
lib_pseudo = {"atom type":{"loc":"specified basis set library","key":"basis set
for atom type"},"atom type2":...}
```
- Each line element is separated by a space
- Write the filenames and paths to the input files next to the corresponding variable (must be named sew0\_file, psd\_file e.t.c)
- Give the fragment library and pseudopotential library in the [python dictionary format](https://www.w3schools.com/python/python_dictionaries.asp).

Here is an example:

```
filename = example.sew.in
title = example
sew0_file = GdMn2O5_J1.env.sew0
psd_file = GdMn2O5_J1.env.psd
lib_frag = {’Mn’:{’loc’:’’,’key’:’Mn.ano-rcc.Roos.21s15p10d6f4g2h.6s4p3d1f0g.’},
’O’:{’loc’:’’,’key’:’O.ano-rcc.Roos.14s9p4d3f2g.4s3p1d0f’}}
lib_pseudo = {’Gd1’:{’loc’:’PSEUDO’,’key’:’Gd.ECP.Marie.0s.0s.0e-Gd1-GdMn2O5.’},
’Gd2’:{’loc’:’PSEUDO’,’key’:’Gd.ECP.Marie.0s.0s.0e-Gd2-GdMn2O5.’},
’Mn1’:{’loc’:’PSEUDO’,’key’:’Mn.ECP.Marie.0s.0s.0e-Mn1-GdMn2O5.’},
’Mn2’:{’loc’:’PSEUDO’,’key’:’Mn.ECP.Marie.0s.0s.0e-Mn2-GdMn2O5.’},
’O1’:{’loc’:’PSEUDO’,’key’:’O.ECP.Marie.0s.0s.0e-O1-GdMn2O5.’},
’O2’:{’loc’:’PSEUDO’,’key’:’O.ECP.Marie.0s.0s.0e-O2-GdMn2O5.’},
’O3’:{’loc’:’PSEUDO’,’key’:’O.ECP.Marie.0s.0s.0e-O3-GdMn2O5.’},
’O4’:{’loc’:’PSEUDO’,’key’:’O.ECP.Marie.0s.0s.0e-O4-GdMn2O5.’}}
```

## Output

If the setup has been performed correctly, the script will output ‘File has been created’ in the terminal and the sew.in file
will appear in the input file directory.
