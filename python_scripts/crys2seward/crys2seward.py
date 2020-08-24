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


def cell_list(disp_input_file):
    """Generates cell grid from input file in directory priority of:

    1) Quoted relative or full path for cell file
    2) Quoted relative or full path for crystal output file
    3) Current working directory

    From this directory, a parsed list of all cell file names is generated

    Parameters
    ----------
    disp_input_file: str
        Relative, absolute or full path of the disp_solve input file
        containing crystal and cell file outputs as well as designated
        electric field ranges and data input method (direct, auto, charge)

    Returns
    ----------
    grid_dir: str
        Directory of cell file grid
    cell_files: list
        List containing strings detailing the filenames of each cell file
        in the grid
    """

    # This line generates the grid using disp_solve.py and assigns
    # the grid's target directory to the grid_dir variable
    grid_dir = dsolve.read_input(disp_input_file)

    # Cell file list, make sure to only take files and not directories
    # for the list
    cell_files = [f for f in os.listdir(grid_dir)
                  if os.path.isfile(os.path.join(grid_dir, f))]

    return grid_dir, cell_files


def sim_cells(envin, envout, disp_input_file, sew0_file, psd_file):
    """Generates cell grid for cell_list, and runs xenv15 using
    each cell file from the disp_solve.py output as an input.
    Each output sew0 and psd file is renamed according to the
    displacement field.

    Parameters
    ----------
    envin: str
        Absolute path of envin or env.in file for xenv15
    envout: str
        Absolute path of env.out file for xenv15
    disp_input_file: str
        Full, relative or absolute path of the input file for
        disp_solve.py
    sew0_file: str
        The absolute path of the xenv sew0 output file:
        e.g ymno3_d1.sew0 or TbMn2O5_J2.env.sew0
    psd_file: str
        The absolute path of the xenv psd output file:
        e.g ymno3_d1.psd or TbMn2O5_J2.env.psd
        
    Returns
    ----------
    """

    cell_data = cell_list(disp_input_file)
    grid_dir = cell_data[0]
    cell_list = cell_data[1]

    for i in range(len(cell_list)):
        fname = cell_list[i]
        cell_dir = grid_dir + "/" + fname
        change_cell(envin, cell_dir)
        run_xenv(envin, envout)
        
        field = fname.split(".")[0]
        sew0_name = sew0_file.split(".")
        psd_name = psd_file.split(".")
        
        if sew0_name[-2] != "env":
            sew0_name.insert(-1,"env")

        if psd_name[-2] != "env":
            psd_name.insert(-1,"env")

        sew0_name.insert(-2,field)
        sew0_name = ".".join(sew0_name)

        psd_name.insert(-2,field)
        psd_name = ".".join(psd_name)
        
        os.rename(sew0_file, sew0_name)
        os.rename(psd_file, psd_name)
        print(sew0_name, psd_name)



# Test inputs for bug checking functions
envin = "ymno3_d1.envin"
envout = "ymno3_d1.env.out"
cell_file = "../ymno3.new.cell"
disp_input = "../disp_input"
sew0_file = "ymno3_d1.sew0"
psd_file = "ymno3_d1.psd"

# run_xenv(envin, envout)
# change_cell(envin, cell_file)
# cell_list(disp_input)
sim_cells(envin, envout, disp_input, sew0_file, psd_file)
