import numpy as np
import os
import sys
# env2seward full or relative path goes in sys.path.insert
sys.path.insert(1, '../../../python_scripts/env2seward')
import env2seward as e2s
# Same for disp_solve, generating the scripts as packages
# (requires __init__.py in both directories)
sys.path.insert(1, '../../../python_scripts/disp_solve')
import disp_solve as dsolve


def run_xenv(envin, envout):
    """This function runs xenv15 to generate sew0 and psd output
    files from envin input.

    Parameters
    ----------
    envin: str
        Absolute path of prefix.envin file, since xenv15 only outputs
        files to current working directory
    envout: str
        Absolute path of prefix.env.out file in same cwd

    Returns
    ----------
    None
        Writes results to file
    """

    # Generating command line input
    cmd = "./xenv15 <%s> %s" % (envin, envout)

    # Writing command to terminal
    os.system(cmd)


def change_cell(envin, cell_file):
    """Changes cell file path within envin or env.in files

    Parameters
    ----------
    envin: str
        Path to envin or env.in file that will be edited
    cell_file: str
        Path to cell_file to be written into the envin file

    Returns
    ----------
    None
        Edits existing file
    """

    # Open the envin file for reading
    with open(envin, 'r') as file:
        env = file.readlines()

    # Editing the cell file absolute path
    for i in range(len(env)):

        # Empty line = ignore
        if not len(env[i].split()):
            pass

        # Covering both formats
        elif env[i].split('=')[0] == "NFI_In":
            env_split = env[i].split('"')
            env_split[1] = cell_file
            env[i] = '"'.join(env_split) + '\n'
            ch_index = i
            break

        elif env[i].split()[0] == "fcell":
            env[i+1] = cell_file + '\n'
            ch_index = i+1
            break

        elif i == len(env)-1:
            sys.exit("cell file variable not found")

    with open(envin, 'w') as file:
        file.writelines(env)

# Test inputs for bug checking functions
envin = "ymno3_d1.envin"
envout = "ymno3_d1.env.out"
cell_file = "../ymno3.new.cell"

# run_xenv(envin, envout)
# change_cell(envin, cell_file)
