import os
import sys
import shutil
import numpy as np
# env2seward full or relative path goes in sys.path.insert
sys.path.insert(1, '../../../python_scripts/env2seward')
import env2seward as e2s

# Same for disp_solve, generating the scripts as packages
# (requires __init__.py in both directories)
sys.path.insert(1, '../../../python_scripts/disp_solve')
import disp_solve as dsolve

from pathlib import Path


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


def change_cell(envin, cell_file, cell_init):
    """Changes NFI_In or fcell (cell file path) value within
    envin or env.in files as well as nch (number of atoms in
    the unit cell), posmag (magnetic atom coordinates) and
    potentially atom0 if Marie or Elisa add to this function.

    Parameters
    ----------
    envin: str
        Path to envin or env.in file that will be edited
    cell_file: str
        Path to new (assumedly displaced) cell_file to be
        written into the envin file
    cell_init: str
        Path to initial cell_file to track original magnetic
        atom positions

    Returns
    ----------
    new_envin: str
        Name of newly formed envin file with edited variables
    """

    # Open the envin, initial cell file and
    # new cell file for reading
    with open(envin, 'r') as file:
        env = file.readlines()

    with open(cell_init, 'r') as file:
        cell = file.readlines()
        no_atoms = len(cell)

    with open(cell_file, 'r') as file:
        new_cell = file.readlines()

    cell_coords = {}
    cell_new_coords = {}

    for i in range(len(cell)):
        cell_split = cell[i].split()
        temp_coords = []

        cell_coords["%d" % i] = {}

        for j in range(1, len(cell_split) - 1):
            temp_coords.append(float(cell_split[j].replace("D", "e")))

        cell_coords["%d" % i]['coords'] = temp_coords

    for i in range(len(new_cell)):
        cell_split = new_cell[i].split()
        temp_coords = []

        cell_new_coords["%d" % i] = {}

        for j in range(1, len(cell_split) - 1):
            temp_coords.append(float(cell_split[j].replace("D", "e")))

        cell_new_coords["%d" % i]['coords'] = temp_coords

    # Generating relative path to try to avoid hitting
    # the Env15 character read limit with a long full path
    rel_cell = os.path.relpath(cell_file, os.getcwd())

    posmag_index = 0
    # Editing the cell file absolute path
    for i in range(len(env)):

        # Empty line = ignore
        if not len(env[i].split()):
            pass

        # Covering both cell file entry formats
        elif env[i].split('=')[0] == "NFI_In":
            env_split = env[i].split('"')
            env_split[1] = rel_cell
            env[i] = '"'.join(env_split)
            ch_index = i

        # Second format
        elif env[i].split()[0] == "fcell":
            env[i+1] = cell_file + '\n'
            ch_index = i+1

        # Number of unit cell atoms (format 1)
        elif env[i].split()[0] == "nch":
            env[i+1] = " %d\n" % no_atoms

        # Number of unit cell atoms (format 2)
        elif env[i][:4] == "nch=":
            env[i] = "nch=%d\n" % no_atoms

        # Here you could include a pos0 coordinate-editing if
        # condition, but I did not have enough time for this

        # Editing magnetic atom position displacements, matching them
        # to coordinates in the initial unedited cell file
        if env[i][:6] == "posmag":
            posmag_index += 1

            # Generating coordinates read from envin posmag line
            temp_coords = []
            coords_split = env[i].split("=")[1:][0].split()[0:3]
            for j in range(len(coords_split)):
                temp_coords.append(float("".join([k.replace("d", "e")
                                   for k in coords_split[j] if k != ","])))

            # print("Posmag %d coords:" % posmag_index)
            # print(temp_coords)

            # Reading the cell file coordinates and comparing their
            # separation. Add associated displacement if sep < 1e-1
            for index in cell_coords:
                dict_coords = cell_coords[index]['coords']
                separation_squared = 0

                # print("Initial cell file coordinates")
                # print(dict_coords)

                for j in range(len(dict_coords)):
                    separation_squared += (dict_coords[j]
                                           - temp_coords[j])**2

                separation = np.sqrt(separation_squared)
                # print("Separation: ", separation)

                if separation <= 1e-1:
                    # print("Coordinate match found for posmag %d"
                    #       % (posmag_index + 1))
                    match_index = int(index)
                    replace_coords = cell_new_coords[index]['coords']
                    break

            # Generating new coordinates for the new envin file
            replace_formatted = ", ".join(['{0:.12f}'.format(i)+"d0"
                                           for i in replace_coords])

            replace_line = "=".join([env[i].split("=")[0],
                                    replace_formatted])+"\n"

            for j in range(len(replace_line)):
                if replace_line[j] == "=":
                    replace_line = replace_line[:j+1] + " "\
                                   + replace_line[j+1:]
                    break

            env[i] = replace_line
            # print("new posmag line: ")
            # print(env[i])

    # Generating new envin file name and writing the new envin file
    envin_split = envin.split(".")
    envin_split.insert(-1, "new")
    new_envin = ".".join(envin_split)

    with open(new_envin, 'w') as file:
        file.writelines(env)

    return new_envin


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

    cell_files.sort()

    return grid_dir, cell_files


def sim_cells(envin, envout, disp_input_file, sew0_file,
              psd_file, cell_init):
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
    cell_init: str
        Path of initial cell file to be used as input for xenv

    Returns
    ----------
    directory: str
        Full path of folder containing the sew0/psd grid
    """

    # Defining old grid directory and list of cell files

    print("###################Generating grid#############################")

    cell_data = cell_list(disp_input_file)

    grid_dir = cell_data[0]
    list_cell = cell_data[1]

    # Name of new grid folders for sew0 and psd files
    grid_sew0 = grid_dir.split("/")[-1].replace("_cell", "_sew0")
    directory_sew0 = os.getcwd() + "/" + grid_sew0

    grid_psd = grid_dir.split("/")[-1].replace("_cell", "_psd")
    directory_psd = os.getcwd() + "/" + grid_psd

    # Generating new grid folder
    Path(directory_sew0).mkdir(parents=True, exist_ok=True)
    Path(directory_psd).mkdir(parents=True, exist_ok=True)

    # Running xenv15 after changing the envin file according to
    # the new displaced unit cell
    for i in range(len(list_cell)):
        fname = list_cell[i]
        cell_dir = grid_dir + "/" + fname
        new_envin = change_cell(envin, cell_dir, cell_init)
        run_xenv(new_envin, envout)

        field = ".".join(fname.split(".")[1:-1])

        # Generating new names for the output files
        sew0_name = sew0_file.split(".")
        psd_name = psd_file.split(".")

        if sew0_name[-2] != "env":
            sew0_name.insert(-1, "env")

        if psd_name[-2] != "env":
            psd_name.insert(-1, "env")

        sew0_name.insert(-2, field)
        sew0_name = ".".join(sew0_name)

        psd_name.insert(-2, field)
        psd_name = ".".join(psd_name)

        # Renaming output files
        os.rename(sew0_file, sew0_name)
        os.rename(psd_file, psd_name)

        print("#####################################################")
        print("xenv15 generated files: %s, %s" % (sew0_name, psd_name))

        # Moving output files to the desired folder
        shutil.move(os.path.join(os.getcwd(), sew0_name),
                    os.path.join(directory_sew0, sew0_name))
        shutil.move(os.path.join(os.getcwd(), psd_name),
                    os.path.join(directory_psd, psd_name))

    return directory_sew0, directory_psd


def sew_in_grid(envin, envout, disp_input_file, env2sew_input_file,
                sew0_file, psd_file, cell_init):
    """This function uses the grid of xenv15 outputs (prefix.env.sew0
    and prefix.env.psd) to generate a grid of corresponding
    prefix.sew.in files

    Parameters
    ----------
    envin: str
        Absolute path of envin or env.in file for xenv15
    envout: str
        Absolute path of env.out file for xenv15
    disp_input_file: str
        Full, relative or absolute path of the input file for
        disp_solve.py
    env2sew_input_file: str
        Full, relative or absolute path of the input file for
        env2seward.py
    sew0_file: str
        The absolute path of the xenv sew0 output file:
        e.g ymno3_d1.sew0 or TbMn2O5_J2.env.sew0
    psd_file: str
        The absolute path of the xenv psd output file:
        e.g ymno3_d1.psd or TbMn2O5_J2.env.psd
    cell_init: str
        Path of initial cell file to be used as input for xenv

    Returns
    ----------
    None
        sew.in file grid
    """

    # Directory of sew0 and psd grids
    new_grid_dirs = sim_cells(envin, envout, disp_input_file, sew0_file,
                              psd_file, cell_init)

    dir_sew0 = new_grid_dirs[0]
    dir_psd = new_grid_dirs[1]

    for i in range(len(dir_sew0)):
        if dir_sew0[i:i+2] == "Ex":
            start = i

        elif dir_sew0[i:i+4] == "sew0":
            end = i-1

    # Naming and generating new directory for the sew.in grid
    dir_sewin = dir_sew0[start:end]+"_sew.in"

    Path(dir_sewin).mkdir(parents=True, exist_ok=True)

    # List of sew0 and psd files generated in the dir_sew0
    # and dir_psd directories
    sew_list = [f for f in os.listdir(dir_sew0)
                if os.path.isfile(os.path.join(dir_sew0, f))]
    sew_list.sort()

    psd_list = [f for f in os.listdir(dir_psd)
                if os.path.isfile(os.path.join(dir_psd, f))]
    psd_list.sort()

    # Reading env2seward input for editing
    with open(env2sew_input_file, 'r') as file:
        e2s_input = file.readlines()

    # Making a new env2seward input file to process each sew0 and psd file
    for i in range(len(sew_list)):
        sewin_name = sew_list[i].replace(".env.sew0", ".sew.in")

        for j in range(len(e2s_input)):
            e2s_split = e2s_input[j].split()

            if e2s_split[0] == "filename":
                e2s_split[-1] = sew_list[i].split('.')[0] + ".sew.in"
                e2s_input[j] = " ".join(e2s_split) + "\n"
                sewin_file = e2s_split[-1]

            elif e2s_split[0] == "sew0_file":
                e2s_split[-1] = dir_sew0 + "/" + sew_list[i]
                e2s_input[j] = " ".join(e2s_split) + "\n"

            elif e2s_split[0] == "psd_file":
                e2s_split[-1] = dir_psd + "/" + psd_list[i]
                e2s_input[j] = " ".join(e2s_split) + "\n"

        with open('new_e2s_input', 'w') as file:
            file.writelines(e2s_input)

        print("#####################################################")
        print("env2seward generated output file: %s" % sewin_name)

        # Writing new sew.in file from new env2seward input file
        e2s.fileinput('new_e2s_input')

        # Renaming the sew.in file
        os.rename(sewin_file, sewin_name)

        # Moving it to the desired directory
        shutil.move(os.path.join(os.getcwd(), sewin_name),
                    os.path.join(os.path.join(os.getcwd(),
                                 dir_sewin), sewin_name))

    print("#####################################################")
    print("sew.in grid generated")


def read_input_file(c2s_input):
    """This function reads the master input file containing info
    required for sew_in_grid, disp_solve and env2seward inputs

    Parameters
    ----------
    c2s_input: str
        Full, relative or absolute path to input file containing
        paths to envin, envout, initial cell file, disp_solve input
        file, names of xenv15 system output files and the contents
        of the disp_solve and env2seward input files, separated for
        readability

    Returns
    ----------
    None
        Outputs sew.in grid to the current working directory
    """

    # Opening c2s input file for reading
    with open(c2s_input, 'r') as file:
        c2s = file.readlines()

    # Parse c2s input file for variables
    for i in range(len(c2s)):
        c2s_split = c2s[i].split()

        if c2s_split[0] == "envin":
            envin = " ".join(c2s_split[2:])

        if c2s_split[0] == "envout":
            envout = " ".join(c2s_split[2:])

        if c2s_split[0] == "cell_init":
            cell_init = " ".join(c2s_split[2:])

        if c2s_split[0] == "xenv_sew0":
            xenv_sew0 = " ".join(c2s_split[2:])

        if c2s_split[0] == "xenv_psd":
            xenv_psd = " ".join(c2s_split[2:])

    # Calling grid generation function with the parsed variables
    sew_in_grid(envin, envout, c2s_input, c2s_input,
                xenv_sew0, xenv_psd, cell_init)


def command_line():
    """Allows for command line input - pass input path as:
    python3 <crys2seward.py_path> <input_file_path>

    No parameters or returns.
    """
    print("CONSIDER THE FOLLOWING QUESTIONS CAREFULLY:")
    basis_bool = e2s.ask_user("Is the system's crystallographic basis"+\
                              " and primitive basis the same?")
    atom_order_bool = e2s.ask_user("Are the displacements in the same"+\
                                   " order as the cell file atoms? Is"+\
                                   " there a one to one correspondence?")
    if len(sys.argv) >= 2 and basis_bool and atom_order_bool:
        read_input_file(sys.argv[1])
    
    else:
        sys.exit("Do not run the code if your answers to these questions"+\
                 " are no.")

# Test inputs for bug checking functions

# envin = "ymno3_d1.envin"
# envout = "ymno3_d1.env.out"
# cell_file = "../ymno3.new.cell"
# disp_input = "../disp_input"
# sew0_file = "ymno3_d1.sew0"
# psd_file = "ymno3_d1.psd"
# e2s_input = "ymno3_d1.new.in"
# c2s_input = "../ymno3_d1.c2s.in"
# envin2 = "CuO_P1.envin"
# cell_file2 = "CuO_P1.loto.cell"
# cell_new = "Ex_0.0_0.15_0.15__Ey_0.0_0.15_0.15__Ez_0.0_0.15_0.15_cell"+\
#            "/CuO.loto.0.15_0.15_0.15.cell"

# run_xenv(envin, envout)
# change_cell(envin2, cell_new, cell_file2)
# cell_list(disp_input)
# sim_cells(envin, envout, disp_input, sew0_file, psd_file, cell_file)
# sew_in_grid(envin, envout, c2s_input, c2s_input,
#             sew0_file, psd_file, cell_file)
# read_input_file(c2s_input)

command_line()
