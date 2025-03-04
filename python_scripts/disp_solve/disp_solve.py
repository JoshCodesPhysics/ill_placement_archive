#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 08:35:02 2020

@author: joshhorswill10
"""

import numpy as np
import sys
import itertools as it
import copy
import os.path
from tabulate import tabulate
import ast
from pathlib import Path
import json

np.set_printoptions(precision=15)


def num_atoms(born_file):
    """Determines number of atoms in the unit cell from parsing 
    the CRYSTAL output file
    
    Parameters
    ----------
    born_file: str
        The path of the file containing the Born charge data
    
    Returns
    ----------
    atom_num: int
        The number of atoms in the unit cell for the system
    """
    
    born = open(born_file).readlines()
    
    for i in range(len(born)):
        born_split = born[i].split()
        
        if born_split[:5] == ['ATOMS', 'IN', 'THE', 'ASYMMETRIC', 'UNIT']:
            atom_num = int(born_split[-1])
            break
    
    return atom_num


def matrix_data(born_file):
    """Stores Born charge matrix data for each atom as 3x3 matrix corresponding
    to the atom key
    
    Parameters
    ----------
    born_file: str
        The path of the file containing the Born charge data
    
    Returns
    ----------
    matrix_dict: dict
        Dictionary containing 3x3 born charge matrices assigned to
        each atom in the unit cell
    """
    
    atom_num = num_atoms(born_file)
    # Generating atom data and empty storage matrix
    atom_array = np.arange(1, atom_num+1, 1)
    matrix_dict = {"atom%d"%i : np.zeros((3, 3)) for i in atom_array}
    born = open(born_file).readlines()
    
    # Searching for title and truncating relevant data
    for i in range(len(born)):
        born_split = born[i].split()
        if born_split[:3] == ['BORN', 'CHARGE', 'TENSOR.']:
            start = i+1
            end = (i+1) + (7*atom_num)
            break
        
        elif born_split[:4] == ['ATOMIC', 'BORN', 'CHARGE', 'TENSOR']:
            start = i + 2
            end = (i + 2) + (7*atom_num)
            break

    tensor_data = born[start:end]
    
    # Appending data to dictionary as 3x3 matrices
    for i in range(len(tensor_data)):
        tensorsplit = tensor_data[i].split()
        if len(tensorsplit) > 0:
            if tensorsplit[0] == 'ATOM':
                atom_index = int(tensorsplit[1])
                m_index = i+3
                for j in range(3):
                    if j != 0:
                        m_index += 1
                    for k in range(3):
                        if abs(float(tensor_data[m_index].split()[k+1]))\
                                < 1e-12:
                            matrix_dict["atom%d"%atom_index][j, k] = 0e0
                        else:
                            matrix_dict["atom%d"%atom_index][j, k] \
                            = float(tensor_data[m_index].split()[k+1])
    return matrix_dict


def born_tensor(born_file):
    """Appends Born matrix data to a larger matrix where each
    3x3 matrix is a diagonal element of the larger matrix
    
    Parameters
    ----------
    born_file: str
        The path of the file containing the Born charge data
    
    Returns
    ----------
    atom_num: int
        The number of atoms in the unit cell for the system
    """
    
    # Defining matrix data
    m_dict = matrix_data(born_file)
    m_dim = 3*num_atoms(born_file)
    tensor = np.zeros((m_dim, m_dim))
    start_column = 0
    atom_num = 1

    # Looping through tensor columns and adding
    # stepped matrices along the diagonal
    for i in range(m_dim):
        if i != 0 and i % 3 == 0:
            start_column += 3
            atom_num += 1
        
        for j in range(3):
            tensor[i, j+start_column] = float(
                   m_dict['atom%d'%atom_num][i%3, j])
    
    return tensor


def born_dat(born_file):
    """Reads BORN.DAT file from crystal output directly, initally
    written to compare the data between the old displacements and
    new displacements. Outputs (3N, 3N) dimension diagonalised Born
    tensor with reduced error"

    Parameters
    ----------
    born_file: str
        BORN.DAT path to read born tensor data directly

    Returns
    ----------
    born_tensor: array
        Born charge tensor array
    """

    with open(born_file, 'r') as file:
        born_data = file.readlines()

    atom_num = int(len(born_data) / 3)
    atom_array = np.arange(1, atom_num+1, 1)
    born_tensor = np.zeros((3*atom_num,3*atom_num))

    m_number = 0
    for i in range(len(born_data)):
        if i % 3 == 0:
            for j in range(3):
                for k in range(3):
                    if abs(float(born_data[i + j].split()[k])) <= 1e-12:
                        pass
                    
                    else:
                        born_tensor[m_number + j, m_number + k]\
                            = float(born_data[i + j].split()[k])
            
            m_number += 3

    return born_tensor

###################################################################
def parse_hessian(hess_file, thresh):
    """Parses Hessian to find index of pairs and diagonal elements
    
    Parameters
    ----------
    hess_file: str
        The path of the file containing the symmetric Hessian data
    thresh: float
        Maximum difference between two values for a pair to be mapped
    
    Returns
    ----------
    hess_pairs: list
        Paired values in tuples
    hess_diag: list
        Diagonal values with no matches
    hess_copy: list
        All values, unpaired and unsorted with zeros included
    """
    
    hess_data = []
    hess = open(hess_file).readlines()
    
    for i in range(len(hess)):
        hess_split = hess[i].split()
        # Appending all values to hess_data from input file, unless the value
        # is very small, then it is appended as zero
        for j in range(len(hess_split)):
            if abs(float(hess_split[j]))<= 1e-15:
                hess_data.append(0.0)
            
            else:
                # Appending data and keeping the original index for later
                hess_data.append((float(hess_split[j]), i*len(hess_split) + j))
    
    # Keeping a copy of the data with zeros included
    hess_copy = copy.deepcopy(hess_data)
    # Removing zeros
    hess_data = [i for i in hess_data if i != 0.0]
    hess_pairs = []
    hess_diag = []
    flag = True
    
    # While there are still elements in hess_data, append and remove pairs
    # with their respective indexes
    # as well as unique elements with their indexes.
    while flag:
        # Uncomment below statement if you want a 'progress bar'
        # print("initial element length ",len(hess_data),
        # "twice pair length + diag length ",
        # 2*len(hess_pairs)+len(hess_diag))
        if len(hess_data) == 1:
            # Removes bug for last element
            hess_diag.append(hess_data[0])
            hess_data.remove(hess_data[0])
            break

        # Break once list is empty
        if len(hess_data) == 0:
            flag = False
            break
        
        for i in range(1, len(hess_data)):
            # If the difference between two elements satisfies the threshold
            # Add them and their indexes to a pair list and remove both 
            # from hess_data
            if abs(hess_data[0][0]-hess_data[i][0]) <= thresh:
                mean = (hess_data[0][0]+hess_data[i][0])/2.0
                hess_pairs.append([mean,hess_data[0][1], hess_data[i][1]])
                rem1 = hess_data[0]
                rem2 = hess_data[i]
                hess_data.remove(rem1), hess_data.remove(rem2)
                break
            
            # If there are no pairs for this value, add it to diagonal list
            # with it's index and remove
            elif i == len(hess_data)-1:
                hess_diag.append([hess_data[0][0], hess_data[0][1]])
                hess_data.remove(hess_data[0])
    
    # return triple of pairs, diagonal and original data lists
    return hess_pairs, hess_diag, hess_copy


def print_indexes(hess_file, thresh):
    """Prints all paired and unique values with their indexes

    Parameters
    ----------
    hess_file: str
       Path for file containing symmatrix Hessian data
    thresh: float
       Threshold for maximum difference between pairs for them
       to be coupled

    Returns
    ----------
    None
    """
    
    print("Pair indexes")
    print([i[1:] for i in parse_hessian(hess_file, thresh)[0]])
    print("Diag indexes")
    print([i[1] for i in parse_hessian(hess_file, thresh)[1]])


def print_diag(hess_file, thresh, rev):
    """Writes diagonal elements to a text file in either ascending
    or descending order depending on rev boolean
    
    Parameters
    ----------
    hess_file: str
       Path for file containing symmatrix Hessian data
    thresh: float
       Threshold for maximum difference between pairs for them
       to be coupled
    rev: bool
       True = descending diagonal elements in file
       False = ascending

    Returns
    ----------
    None

    """
    
    with open("diagonal.txt",'w') as diag:
        sort_diag = parse_hessian(hess_file, thresh)[1]
        sort_diag = [i[0] for i in sort_diag]
        # Sorting ascending or descending, rev = True creates descending
        sort_diag.sort(reverse=rev)
        print("Length ", len(sort_diag))
        for item in sort_diag:
            diag.write("%s\n"%str(item))


def num_diag(hess_file, nthresh, threshstart):
    """Prints a dictionary containing number of diagonal elements for a range 
    of thresholds nthresh gives number of thresholds to test,
    and threshstart gives starting value that gets smaller and smaller.
    
    Parameters
    ----------
    hess_file: str
       Path for file containing symmatrix Hessian data
    nthresh: int
        Number of thresholds to measure diagonal element number
    threshstart: float
       Starting value for threshold to be made increasingly smaller
       for every input into parse_hessian

    Returns
    ----------
    None

    """
    
    thresh_range = np.ones((nthresh))
    # Generating threshold values
    for i in range(len(thresh_range)):
        if i == 0:
            thresh_range[i] = threshstart
        
        else:
            # One order of magnitude between each index
            thresh_range[i] = 1e-1*thresh_range[i-1]
    
    # Empty dictionary for writing
    diag_dict = {str(i) : None for i in thresh_range}
    
    # Loop writing all diagonal list lengths for given thresholds
    for i in range(len(thresh_range)):
        diag_dict['%s'%str(thresh_range[i])] =\
                  len(parse_hessian("Pm.rePhonons.2123667.HESSFREQ.DAT",
                  thresh_range[i])[1])

    print("Format is 'threshold':number of diagonal elements")
    print(diag_dict)
################################################################

def convert_coords(output_file):
    """Finds conversion matrix and transposes it. Also parses for
    unit cell fractional coordinate data
    
    Parameters
    ----------
    output_file: str
        Path for CRYSTAL output file containing conversion matrix
        and unit cell fractional coordinate data

    Returns
    ----------
    conv_matrix: array
        Conversion matrix parsed from CRYSTAL ouput file and transposed
    frac_data: list
        Fractional coordinate unit cell data parsed from CRYSTAL output
        separated by line
    lat_param_data: str
        Line containing all lattice parameter data

    """

    # Allocating number of unit cell atoms and conversion matrices
    atom_num = num_atoms(output_file)
    out = open(output_file).readlines()
    conv_matrix = np.zeros((3, 3))
    conv_inverse = np.zeros((3, 3))

    # Angstrom to Bohr conversion factor (if needed)
    # ANG2AT = float(1.889725989) 

    # Parsing for data strips
    for i in range(len(out)):
        out_split = out[i].split()

        if out_split == ['DIRECT', 'LATTICE', 'VECTORS', 'CARTESIAN',
                         'COMPONENTS', '(ANGSTROM)']:
            start = i+2
            end = start+3

        elif out_split == ['COORDINATES', 'OF', 'THE', 'EQUIVALENT',
                           'ATOMS', '(FRACTIONARY', 'UNITS)']\
             or out_split == ['COORDINATES', 'OF', 'THE', 'EQUIVALENT',
                              'ATOMS', '(FRACTIONAL', 'UNITS)']:
            start2 =i+4

        elif out_split[:5] == ['NUMBER', 'OF', 'SYMMETRY', 'OPERATORS', ':']:
            end2 = i-1

        elif out_split == ['LATTICE', 'PARAMETERS', '(ANGSTROMS', 'AND',
                           'DEGREES)', '-', 'PRIMITIVE', 'CELL']:
            start3 = i+2

        # elif out_split == ['*', 'ATOM', 'X(ANGSTROM)', 'Y(ANGSTROM)',
        #                     'Z(ANGSTROM)']:
        #    start3 = i+2
        #    end3 = start3+atom_num
        #    break

    # Conversion matrix strip
    conv_data = out[start:end]

    # Unit cell in fractional units strip
    if "start2" not in locals():
        for i in range(len(out)):
            out_split = out[i].split()

            if out_split[:5] == ['ATOMS', 'IN', 'THE', 'ASYMMETRIC', 'UNIT']:
                start2 = i+3
                end2 = start2 + atom_num
                break

        frac_data = out[start2 : end2]
    else:
        frac_data = out[start2 : end2]
    
    # Cartesian coordinates to test conversion matrix
    # test_data = out[start3:end3]

    # Finding lattice parameter data
    lat_param_data = out[start3]

    # Appending conversion matrix
    for i in range(len(conv_data)):
        conv_split = conv_data[i].split()
        for j in range(len(conv_split)):
            conv_matrix[i,j] = float(conv_split[j])

    # Transposing it so we can solve for u,v,w
    conv_matrix = np.transpose(conv_matrix)
    # conv_inverse = np.linalg.inv(conv_matrix)

    return [conv_matrix, frac_data, lat_param_data]


def conv_hessian(hess_matrix, output_file):
    """Converts Hessian from fractional units to atomic units

    Parameters
    ----------
    hess_matrix: array
        Contains non-converted Hessian matrix, symmetrised
        along the diagonal
    output_file: str
        Contains path to CRYSTAL output file containing
        necessary conversion matrix and lattice parameter data
    
    Returns:
    ----------
    hess_matrix: array
        Contains Hessian matrix converted to fractional atomic units?
        In the a,b,c lattice vector basis? Yes this may be changed.
        Not sure yet.
    """
 
    ANG2BOHR = float(1.889725989) 
    conv_data = convert_coords(output_file)
    conv_matrix = conv_data[0]
    
    # Convert lattice parameters into atomic units for
    # matrix element division
    lat_param = [float(ANG2BOHR*float(i))
                 for i in conv_data[2].split()[:3]]
    mdim = len(hess_matrix)
    
    for i in range(mdim):
        for j in range(mdim):
            # if i % 3 == 0:
            #    hess_matrix[i : i+3, j] = np.matmul(conv_matrix,
            #    hess_matrix[i : i+3, j])
            # if j % 3 == 0:
            #    hess_matrix[i, j : j+3] = np.matmul(conv_matrix,
            #    hess_matrix[i, j : j+3])
            # if i % 3 == 0:
            #    if j % 3 == 0:
            #        hess_matrix[i : i+3, j : j+3] = np.matmul(conv_matrix,
            #        hess_matrix[i : i+3, j : j+3])
            hess_matrix[i, j] = (hess_matrix[i, j]/(float(lat_param[i%3])
                                                   *float(lat_param[j%3])))
    return hess_matrix


def hessian(born_file, hess_file):
    """Returns Hessian data from input file as a square symmetric matrix

    Parameters
    ----------
    born_file: str
        Path of file containing born charge data
    hess_file: str
        Path of file containing Hessian data

    Returns
    ----------
    l_matrix: array
        Square L matrix containing values appended from hess_file
        Top right values are all zero
    conv_hess: array
        Processed Hessian, mirrored so it is symmetric. Converted
        from Hartree energy to fractional a.u., still in a,b,c
        basis
    """
    
    # Determines dimensions of Hessian
    m_dim = 3*num_atoms(born_file)
    # Empty square matrix
    tensor_hess = np.zeros((m_dim, m_dim))
    hess = open(hess_file).readlines()
    final_block = False
    start_found = False

    # Finding relevant data range
    for i in range(len(hess)):
        hess_split = hess[i].split()

        # Start condition for Pm.rePhonons file
        if hess_split == ['LOW', 'HALF', 'OF', 'HRED', 'AFTER',
                'SYMMETRISATION']:
            start = i + 4
            start_found = True

        # Start condition for Pm.reOptGeom file
        elif hess_split == [str(i) for i in range(1, len(hess_split)+1)]\
                           and len(hess_split) > 0\
                           and hess[i-1] == '\n'\
                           and len(hess[i+1].split()) == 2\
                           and start_found == False:
            start = i
            start_found = True
        
        # Detecting final block to prepare for end truncation
        if ((str(m_dim) in hess_split\
                and str(m_dim - 1) in hess_split)\
                or (str(m_dim) in hess_split and str(m_dim+1) in hess_split)\
                or (hess_split == [str(m_dim)])):
            if start_found:
                final_block = True

        # End detected, break the loop
        if len(hess_split) > 0:
            if hess_split[0] == str(m_dim)\
                    and final_block\
                    and start_found:
                end = i+1
                break

    # Defining new range
    hess_data = hess[start:end]
    # Generalised column numbers for multiple format possibilities
    ncolumns = len(hess_data[0].split())
    # First column starts at 0
    # Increase by number of columns for every data block
    first_column = 0
    
    # If we are in the last block and we hit the last three rows and columns:
    # break
    last_block = False
    for i in range(len(hess_data)):
        # Splitting lines into individual elements with no spacing
        rowsplit = hess_data[i].split()

        # If we hit an empty line, increase the column number by ncolumns
        if len(rowsplit) == 0:
            first_column += ncolumns
            # print("empty space")

        # If we hit a row of column labels,
        # do nothing except for if it is the last block
        elif rowsplit[:2] == [str(first_column+1), str(first_column+2)]:
            # print("Column label row")
            if str(int(m_dim+1)) in rowsplit:
                last_block = True

        # Data lines
        else:
            # Row number is first element given
            row = int(rowsplit[0]) - 1
            # If we hit last block and last three rows, stop the loop
            if last_block == True and row >= m_dim:
                break
            # If we hit the last three rows in a non-final block,
            # don't take the data
            elif last_block == False and row >= m_dim:
                pass
            # Appending values to tensor
            else:
                for j in range(1, len(rowsplit)):
                    tensor_hess[row, j+first_column-1] = float(rowsplit[j])

    # Copy of L matrix if it is needed
    l_matrix = copy.deepcopy(tensor_hess)

    # Symmetrising
    for i in range(m_dim):
        for j in range(m_dim):
            if i >= j:
                if abs(tensor_hess[i, j]) <= 1e-12:
                    tensor_hess[i, j] = 0
            else:
                if abs(tensor_hess[j, i]) <= 1e-12:
                    tensor_hess[i, j] = 0
                else:
                    tensor_hess[i, j] = float(tensor_hess[j, i])

    # We don't need to conver the Hessian anymore, it is already
    # in the Cartesian basis in a.u.
    # conv_hess = conv_hessian(tensor_hess, hess_file)

    # Add [0] suffix to function for L matrix and [1] for symmetric matrix
    return l_matrix, tensor_hess


def hess_dat(born_file, hess_file, output_file):
    """This function generates the Hessian array from the system
    HESSIAN.DAT file, initially written to compare the old and
    new displacement results

    Parameters
    ----------
    born_file: str
        Path to Born file to extract number of atoms in the
        unit cell (matching dimensions)
    hess_file: str
        Path to HESSIAN.DAT file to extract array

    Returns
    ----------
    conv_hess: array
        Converted Hessian matrix to fractional a.u. in a,b,c
        lattice parameter basis
    """
    with open(born_file, 'r') as file:
        dim = len(file.readlines())

    with open(hess_file, 'r') as file:
        hess_data = file.readlines()

    tensor_hess = np.zeros((dim, dim))
    temp_list = []
    
    for i in range(len(hess_data)):
        hess_split = hess_data[i].split()
        
        for j in range(len(hess_split)):
            temp_list.append(float(hess_split[j]))

    row = 0
    for i in range(len(temp_list)):
        if i % dim == 0 and i != 0:
            row += 1
        tensor_hess[row, i%dim] = temp_list[i]

    for i in range(dim):
        for j in range(dim):
            if i >= j:
                if abs(tensor_hess[i, j]) <= 1e-12:
                    tensor_hess[i, j] = 0

            else:
                if abs(tensor_hess[i, j]) <= 1e-12:
                    tensor_hess[i, j] = 0
                
                else:
                    tensor_hess[i, j] = tensor_hess[j, i]
    # conv_hess = conv_hessian(tensor_hess, output_file)

    return tensor_hess


def born_input(born_file):
    """Checks if born input file is sourced from the crystal output
    such as the .loto.out filetype, or from the .DAT files

    Parameters
    ----------
    born_file: str
        Path of born input file
    
    Returns
    ----------
    b_tensor: array
        Appropriate function used to generated born tensor
    """

    if ".DAT" in born_file:
        b_tensor = born_dat(born_file)

    else:
        b_tensor = born_tensor(born_file)
    
    return b_tensor

def hess_input(born_file, hess_file, output_file):
    """Checks if Hessian input file is sourced from the crystal output
    such as the .loto.out filetype, or from the .DAT files. If Hessian
    input is .DAT, then so must the Born input for this function to
    work

    Parameters
    ----------
    born_file: str
        Path of born input file
    hess_file: str
        Path of Hessian input file
    output_file: str
        Path of crystal output file (e.g .loto files)

    Returns
    ----------
    tensor_hess: array
        Appropriate function used to generate Hessian matrix
    """

    if ".DAT" in hess_file:
        tensor_hess = hess_dat(born_file, hess_file, output_file)

    else:
        tensor_hess = hessian(born_file, hess_file)[1]

    return tensor_hess

########################################################################
def convert_coords2(out_file):
    """Parses for lattice parameters and converts them into unit cell edge
    vectors in cartesian coordinates

    Parameters
    ----------
    out_file: str
        CRYSTAL output file to be parsed for data to supplement non-general
        conversion of a unit cell

    Returns
    ----------
    lattice_vectors: list
        List of lattice vectors generated from lattice parameter data
    unit_vectors: list
        List of unit lattice vectors to be used for a basis
    """
    
    # This function isn't very general so is not used
    out = open(out_file).readlines()
    lat_list = []
    for i in range(len(out)):
        out_split = out[i].split()
        if out_split == ['LATTICE', 'PARAMETERS', '(ANGSTROMS', 'AND',\
                'DEGREES)', '-', 'CONVENTIONAL', 'CELL']:
            lat_data = out[i+2].split()
            break
    
    for i in range(len(lat_data)):
        lat_list.append(float(lat_data[i]))
    
    a = lat_list[0]
    b = lat_list[1]
    c = lat_list[2]
    
    alpha = np.radians(lat_list[3])
    beta = np.radians(lat_list[4])
    gamma = np.radians(lat_list[5])
    
    a_vector = np.array([a, 0.0, 0.0])
    b_vector = np.array([b*np.cos(gamma), b*np.sin(gamma), 0.0])
    
    omega = a*b*c*np.sqrt(1-np.cos(alpha)*np.cos(alpha)-np.cos(beta)
            *np.cos(beta)-np.cos(gamma)*np.cos(gamma)
            +2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
    
    c_vector = np.array([c*np.cos(beta),
               c*((np.cos(alpha)-np.cos(beta)*np.cos(gamma))/(np.sin(gamma))),
               omega/(a*b*np.sin(gamma))])
    
    lattice_vectors = np.array([a_vector, b_vector, c_vector])
    
    unit_vectors = np.array([float(vec/np.sqrt(vec.dot(vec))) 
        for vec in lattice_vectors])
    return lattice_vectors, unit_vectors
###############################################################

def evec(ex, ey, ez, mdim):
    """Generates electric field vector from kV/m a,b,c inputs and
    the converts into a.u., dividing by corresponding 
    lattice parameters as conversion*(Ex,Ey,Ez)
    
    Parameters
    ----------
    ex: float
        Input for electric field in 'x' vector axis (Ex)
    ey: float
        Input for electric field in 'y' vector axis (Ey)
    ez: float
        Input for electric field in 'z' vector axis (Ez)
    mdim: int
        Number of rows in vector / dimension of qE = -Hd equation

    Returns:
        E vector in fractional atomic units with mdim sets of
        Ex, Ey, Ez coordinates (uniform field)
    """
    
    ATOMIC_FIELD = float(5.14220674763e11)
    KV = float(1e3)
    CONVERSION = float((KV/(ATOMIC_FIELD)))
    
    e_vector = np.zeros((mdim))
    
    for i in range(mdim):
        if i % 3 == 0:
            e_vector[i] = float(CONVERSION*ex)
        
        elif i % 3 == 1:
            e_vector[i] = float(CONVERSION*ey)
        
        else:
            e_vector[i] = float(CONVERSION*ez)
    
    return e_vector


def save_hess_born(born_file, hess_file, output_file, ex, ey, ez):
    """Writes arrays to binary files for fortran processing
    Enter electric field components in kV/m
    
    Parameters
    ----------
    born_file: str
        Path of file containing Born charge data
    hess_file: str
        Path of file containing Hessian data
    output_file: str
        Path of crystal output file containing unit cell
        data
    ex: float
        Input for electric field in 'x' vector axis (Ex)
    ey: float
        Input for electric field in 'y' vector axis (Ey)
    ez: float
        Input for electric field in 'z' vector axis (Ez)

    Returns
    ----------
    None
        Writes Hessian, Born and E vector binary files 
        for Fortran script input
    """
 
    h_tensor = hess_input(born_file, hess_file, output_file)
    q_tensor = born_input(born_file)
    e_vector = evec(ex, ey, ez, len(h_tensor))

    # Generating e_vector from inputs
    dim_array = np.array([len(h_tensor)])

    # Saves all matrices and dimension to binary files
    # to be read by lapack_solve.f03
    h_tensor.tofile('hess_matrix')
    q_tensor.tofile('born_matrix')
    e_vector.tofile('e_vector')
    dim_array.tofile('n_array')

    # np.save('hess_matrix',h_tensor,allow_pickle=True,fix_imports=False)
    # np.save('born_matrix',q_tensor,allow_pickle=True,fix_imports=False)
    print("q,E,H matrices written to binary files")


def solve_equation(born_file, hess_file, output_file, ex, ey, ez):
    """Generates all arrays necessary for qE = -Hd equation and
    solves for d (displacement). Equation occurs in a.u. and the
    displacements are in the same units in the same Cartesian basis

    Parameters
    ----------
    born_file: str
        Path of file containing Born charge data
    hess_file: str
        Path of file containing Hessian data
    output_file: str
        Path of crystal output file containing unit cel
        data
    ex: float
        Input for electric field in 'x' vector axis (Ex)
    ey: float
        Input for electric field in 'y' vector axis (Ey)
    ez: float
        Input for electric field in 'z' vector axis (Ez)

    Returns
    ----------
    displacement: array
        Contains displacments of each unit cell atom in a.u.
        Each group of three numbers e.g index 0,1,2 is a set
        of 3D displacement coordinates in a,b,c space
    """

    q_tensor = born_input(born_file)
    h_tensor = hess_input(born_file, hess_file, output_file)
    m_dim = len(h_tensor)

    # Generating e_vector from inputs (uniform field)
    e_vector = evec(ex, ey, ez, m_dim)

    # Generating matrices for linalg solve
    a = h_tensor
    b = np.matmul(q_tensor, e_vector)
    displacement = np.linalg.solve(a, b)

    # save_hess_born(born_file,hess_file,lat_param,ex,ey,ez)
    
    # Printing components of linear equation including the solution
    # print("Hessian")
    # print(a)
    # print("qE product matrix: ")
    # print(b)
    # print("Born tensor: ")
    # print(q_tensor)
    # print("Electric field: ")
    # print(e_vector)
    # print("Displacement solutions: ")
    # print(displacement)
    # print("Error matrix: ")
    # print(abs(np.matmul(a, displacement)-b))
    return displacement


def convert_disp(born_file, hess_file, output_file, ex, ey, ez):
    """Applies conversion matrix to displacement solutions.
    Operation takes place in Angstrom, so all displacement coordinates
    are converted to Angstrom. Output coordinates are fractional units
    in the lattice vector basis (a,b,c).

    Parameters
    ----------
    born_file: str
        Path of born file
    hess_file: str
        Path of Hessian file
    output_file: str
        Path of CRYSTAL output file containing the conversion matrix,
        the lattice parameters for the electric field input into
        solve_equation, and the unit cell fractional coordinates
        for output into functions further down the chain.
    ex: float
        Input for electric field in 'x' vector axis (Ex)
    ey: float
        Input for electric field in 'y' vector axis (Ey)
    ec: float
        Input for electric field in 'z' vector axis (Ez)
    
    Returns
    -----------
    convdisp: array
        List of N sets of fractional unit coordinates converted
        from the solve_equation displacements, where N is the number
        of atoms in the unit cell.
    frac: list
        List containing each line of the unit cell coordinates
        so we can later generate a cell file and add the displacements
        to it if needed.
    """
    
    # Calling all necessary data
    BOHR2ANG = float(0.529177210903)
    conv_data = convert_coords(output_file)
    conv_matrix = conv_data[0]
    frac = conv_data[1]
    lat_param_data = conv_data[2]
    lat_param = [float(i) for i in lat_param_data.split()[:3]]
    
    displacement = solve_equation(born_file, hess_file, output_file,
                                  ex, ey, ez)
    convdisp = []
    
    # Applying conversion matrix to one set of x,y,z coordinates at a time
    for i in range(len(displacement)):
        if i == len(displacement) - 1:
            pass
        
        elif i % 3 == 0:
            # Conversion to Angstrom since conversion matrix
            # is in Angstrom units
            xyz = BOHR2ANG*np.array(displacement[i : i+3])
            newcoords = np.linalg.solve(conv_matrix, xyz)
            
            # Appending converted coordinates to array, divided by
            # the lattice parameters for a, b, c now in the
            # crystallographic basis
            for i in range(len(newcoords)):
                convdisp.append(newcoords[i] / lat_param[i])

    convdisp = np.array(convdisp)
    return convdisp, frac


def unit_cell(born_file, hess_file, output_file, cell_file,
              ex, ey, ez, unit_source, *args, **kwargs):
    """Generates the unit cell data in a dictionary
    There are 3 input options:

    1) Takes 'direct' cell file data - use if initial cell file
    matches spatial group/number of atoms in the CRYSTAL output
    unit cell 
    
    2) 'auto'mates data from crystal output completely - If the
    spatial groups do not match, i.e. the initial cell file has
    more or less atoms than the CRYSTAL unit cell, take charge
    data from the initial cell file but use the CRYSTAL unit cell
    to build a cell file to which we can add our converted 
    displacements
    
    3) Automates from crystal output as with 2) but allows
    for 'charge' input from user, so no data is taken from the
    initial cell file

    Moderated by unit_source input (direct,auto,charge) 
    and chosen according to cell_file type

    Parameters
    ----------
    born_file: str
        Path of born file
    hess_file: str
        Path of Hessian file
    output_file: str
        Path of CRYSTAL output file containing the unit cell data
    cell_file: str
        Path of initial cell file containing a form of unit cell
        data and charge data. Can be generated by CRYSTAL or user.
    ex: float
        Input for electric field in 'x' vector axis (Ex)
    ey: float
        Input for electric field in 'y' vector axis (Ey)
    ez: float
        Input for electric field in 'z' vector axis (Ez)
    unit_source: str
        Three key words that determine the source of the initial
        non-disturbed unit cell fractional coordinates and the
        associated charge with each atom type
    *args: dict
        Charge data dictionary e.g {"Y":3.0,"MN":3.0,"O1":-2.0,
        "O2",-2.0}. This is the argument for when unit_source =
        'charge' is chosen, for manual user input

    Returns
    ------
    atom_dict: dict
        Dictionary containing atoms categorised according to the
        irreducible groups, with their respective
        fractional coordinates and charges.
    convdisp: array
        Converted displacements to add to the dictionary coordinates
        later down the function chain.
    """
    
    # Generating preliminary storage lists and required data
    cell = open(cell_file).readlines()
    atom_num = num_atoms(output_file)
    
    conv_data = convert_disp(born_file, hess_file, output_file,
                             ex, ey, ez)
    convdisp = conv_data[0]
    frac_data = conv_data[1]
    atom_nums = [i for i in range(1, atom_num+1)]
    atom_list = []
    
    # Atom key list based on number of atoms in unit cell
    for i in atom_nums:
        atom_list.append("atom %d"%i) 
    
    # Taking data direct from cell file
    if unit_source == "direct":
        # Groups of irreducible atoms
        groups = []
        atom_types = []
        group_num = 1
        start=0
        
        # Appending irreducible atom data
        for i in range(len(cell)):
            cell_split = cell[i].split()
            if cell_split[0] not in atom_types:
                atom_types.append(cell_split[0])
                if i > 0:
                    end = i
                    groups.append(cell[start:end])
                    start = end
            
            if i == len(cell)-1:
                groups.append(cell[start:])
        
        # Dictionary keys:
        group_list = ["group %d"%i for i in range(1, len(groups)+1)]
        # Obtaining charge data for atom_dict directly from cell file
        # No manipulation necessary
        element_charge = {}
        
        for i in range(len(cell)):
            cell_split = cell[i].split()
            element = cell_split[0].upper()
            charge = cell_split[-1]
            if (element, charge) not in element_charge.items():
                element_charge[element]=charge
        
        # Empty dictionary to store important data
        atom_dict = {key:{} for key in group_list}
        count = 1
        
        # Appending dictionary in categories of group number
        # and subcategories of atom number
        for i in range(len(groups)):
            for j in range(len(groups[i])):
                coords_split = groups[i][j].split()
                
                for k in range(len(coords_split)):
                    # convert double precision to single temporarily
                    # for manipulation
                    coords_split[k] = coords_split[k].replace("D+", "e+")\
                                      .replace("D-", "e-")
                
                atom = coords_split[0].upper()
                atom_dict["group %d"%(i+1)]["atom %d"%(count)] =\
                        {"coords":[float(coords_split[k+1])
                            for k in range(3)],"element":atom,
                            "charge":float(element_charge[atom]\
                            .replace("D+", "e+").replace("D-", "e-"))}
                count += 1
    
    else:
        # Else unit data is taken from crystal output file
        # rather than the cell file
        # Another branching will occur when deciding to automate
        # charge data or obtain it manually
        num_groups = len(frac_data) - atom_num + 1

        # Group keys
        group_nums = [i for i in range(1, num_groups+1)]
        gap_index = []

        # Group data
        group_list = []
        
        for i in group_nums:
            group_list.append("group %d"%i)
        
        # Appending indexes of the gaps in the fractional unit cell data
        # to distinguish irreducible atom groups
        for i in range(len(frac_data)):
            frac_split = frac_data[i].split()
            if len(frac_split) == 0:
                gap_index.append(i)
        
        # Appending groups according to these indexes
        groups = []

        if len(gap_index):
            for i in range(len(gap_index)):
                if i == 0:
                    groups.append(frac_data[:gap_index[i]])
                
                elif i == len(gap_index)-1:
                    groups.append(frac_data[gap_index[i-1]+1:gap_index[i]])
                    groups.append(frac_data[gap_index[i]+1:])
                
                else:
                    groups.append(frac_data[gap_index[i-1]+1:gap_index[i]])
       
        # If there are no gaps i.e. the second of the two possible unit
        # cell formats, take a different approach to groups
        else:
            group_list = []
            no_groups = 0
            group_start = 0
            for i in range(len(frac_data)):
                frac_split = frac_data[i].split()
                atom = frac_split[-4]
                # If current atom does not match previous atom, a new
                # group has been started! Append the old one.
                if i == 0:
                    pass

                elif atom != frac_data[i-1].split()[-4]:
                    group_end = i
                    groups.append(frac_data[group_start : group_end])
                    group_start = group_end
                    no_groups += 1

                # We're at the end. Append the final group.
                elif i == len(frac_data) - 1:
                    group_end = i+1
                    groups.append(frac_data[group_start : group_end])
                    no_groups +=1

            group_list = ["group %d" % (i+1) for i in range(no_groups)]

        # print("groups: \n", groups)
        # print("group_list: \n", group_list)
        
        # Automatic charge data generation from initial cell file,
        # potentially from different space group
        # Different spatial group cell file is reason for using auto and
        # charge
        if unit_source == "auto":
            element_charge = {}
            for i in range(len(cell)):
                cell_split = cell[i].split()
                
                if len(cell_split[0]) == 1:
                    element = cell_split[0].upper()

                elif cell_split[0][1].isdigit():
                    element = cell_split[0][0].upper()
                
                else:
                    element = cell_split[0].upper()
                charge = float(cell_split[-1].replace("D","e"))
                
                if (element,charge) not in element_charge.items():
                    element_charge[element] = charge

            # print("element_charge: \n", element_charge)
        
        # Manual charge data input
        elif unit_source == "charge":
            if len(args)>0:
                element_charge = args[0]
        
        # Terminate if incorrect unit_source argument entered
        else:
            sys.exit("Not a recognised argument, see documentation")
        
        atom_dict = {key:{} for key in group_list}
        parsed_atoms = {}
        dist_atoms = {}
        
        # print("length of gap_index: ", len(gap_index))
        # Generating dictionary and then scanning dictionary, relabelling
        # Atoms of same element but different irreducible group
        # e.g O, O -> O1, O2
        for i in range(len(groups)):
            for j in range(len(groups[i])):
                coords_split = groups[i][j].split()
                
                if len(gap_index):
                    atom = coords_split[4].upper()
                    coord_start = 5

                else:
                    atom = coords_split[3].upper()
                    coord_start = 4
                # Charge data included for auto
                if unit_source == "auto":
                    atom_dict["group %d" % (i+1)]["atom %d"
                              %(int(coords_split[0]))] = {"coords":
                              [float(coords_split[coord_start + k])
                              for k in range(3)], "element":atom,
                              "charge":element_charge[atom]}

                # Charge data not included until element data is 
                # reformulated
                elif unit_source == "charge":
                    atom_dict["group %d"%(i+1)]["atom %d"
                              %(int(coords_split[0]))] = {"coords":
                              [float(coords_split[coord_start + k])
                              for k in range(3)],
                              "element":atom}
        
        # Relabel if repeats
        for i in range(len(groups)):
            for j in range(len(groups[i])):
                coords_split = groups[i][j].split()
                if len(gap_index):
                    atom_index = 4
                    atom = coords_split[4].upper()

                else:
                    atom_index = 3
                    atom = coords_split[3].upper()

                if atom not in parsed_atoms and j == 0:
                        parsed_atoms[atom] = i+1
                
                elif coords_split[atom_index] in parsed_atoms:
                    if j > 0:
                        pass
                    
                    elif atom not in dist_atoms:
                        dist_atoms[atom] = 2
                        prev_dict = atom_dict["group %d"
                                              %parsed_atoms[atom]]
                        
                        for k in prev_dict:
                            # print("prev_dict being written for %s"%k)
                            prev_dict[k]["element"] = "%s1"%atom
                        
                        current_dict = atom_dict["group %d"%(i+1)]

                        for k in current_dict:
                            # print("current_dict being written for %s"%k)
                            current_dict[k]["element"] = "%s2"%atom
                    
                    else:
                        dist_atoms[atom] += 1
                        current_dict = atom_dict["group %d"%(i+1)]
                        for k in current_dict:
                            current_dict[k]["element"] = "%s%d"%\
                            (atom, dist_atoms[atom])
        
        # Input manual charge data after correct dictionary is generated
        if unit_source == "charge" and len(args) > 0:
            for group in atom_dict:
                for atom in atom_dict[group]:
                    atom_dict[group][atom]["charge"] =\
                            element_charge[atom_dict[group][atom]\
                                          ["element"]]
        
        # This is for a function later one where we only require
        # the element labels (O1,O2 e.t.c)
        elif unit_source == "charge" and len(args) == 0:
            for group in atom_dict:
                for atom in atom_dict[group]:
                    atom_dict[group][atom]["charge"] = None
    
    # print("##########atom_dict###############")
    # print(atom_dict)
    return atom_dict, convdisp


def modify_cell(born_file, hess_file, output_file, cell_file, dirname,
                ex, ey, ez, unit_source, *args, **kwargs):
    """This function adds converted displacements to the formulated
    unit cell dictionary and writes the results to a new separate
    cell file
    
    Parameters
    ----------
    born_file: str
        Path of born file
    hess_file: str
        Path of Hessian file
    output_file: str
        Path of CRYSTAL output file containing the unit cell data
    cell_file: str
        Path of initial cell file containing a form of unit cell
        data and charge data. Can be generated by CRYSTAL or user.
    dirname: str
        Name of folder created to store new displaced unit cell files
    ex: float
        Input for electric field in 'x' vector axis (Ex)
    ey: float
        Input for electric field in 'y' vector axis (Ey)
    ez: float
        Input for electric field in 'z' vector axis (Ez)
    unit_source: str
        Three key words that determine the source of the initial
        non-disturbed unit cell fractional coordinates and the
        associated charge with each atom type
    *args: dict
        Charge data dictionary e.g {"Y":3.0,"MN":3.0,"O1":-2.0,
        "O2",-2.0}. This is the argument for when unit_source =
        'charge' is chosen, for manual user input

    Returns
    ----------
    None
        Writes edited coordinates to new file
    """
    
    unit_data = unit_cell(born_file, hess_file, output_file, cell_file,
                          ex, ey, ez, unit_source, *args, **kwargs)
    
    atom_dict = unit_data[0]
    convdisp = unit_data[1]
   
    # New file name including electric field parameters
    fsplit = cell_file.split('.')

    if "new" not in fsplit:
        fsplit.insert(-1, '%s_%s_%s'%(str(round(ex, 3)), str(round(ey, 3)),
                      str(round(ez, 3))))
    
    else:
        fsplit = [i.replace("new",'%s_%s_%s'%(str(round(ex, 3)),
                  str(round(ey, 3)), str(round(ez, 3)))) for i in fsplit]
    
    fnew = '.'.join(fsplit).split("/")[-1]

    # Directory of initial cell file and output file
    cell_dir = "/".join(cell_file.split("/")[:-1])
    out_dir = "/".join(output_file.split("/")[:-1])
    
    # If file input is a full/relative path directory rather than filename
    # generate output files in that directory
    if os.path.isdir(cell_dir):
        cd = cell_dir
    
    elif os.path.isdir(out_dir):
        cd = out_dir
    
    else:
        # Else use the current working directory
        cd = os.getcwd()
    
    # Resultant filename given in full path
    completeName = os.path.join(cd+"/"+dirname+"/"+fnew)
    cell_new = open(completeName,'w+')
    
    # Add displacements to unit cell dictionary
    for group in atom_dict:
        for i in range(int(len(convdisp)/3)):
            if "atom %d"%i in atom_dict[group]:
                for j in range(3):
                    atom_dict[group]["atom %d"%i]["coords"][j] +=\
                              convdisp[3*i+j]
    
    # Writing displaced coordinates to a file
    line_list = []
    for group in atom_dict:
        for atom in atom_dict[group]:
            temp_dict = atom_dict[group][atom]
            writelist = [str(temp_dict["element"]),None,None,None,
                    str("{:.15e}".format(temp_dict["charge"]))+'\n']
            for i in range(3):
                writelist[i+1] = str("{:.15e}".format(temp_dict["coords"][i]))
            writeline = "     ".join(writelist)
            line_list.append(writeline.split())
    
    # Converting back to double precision string
    for i in range(len(line_list)):
        for j in range(1,len(line_list[i])):
            line_list[i][j] = line_list[i][j].replace("e+","D+")
            line_list[i][j] = line_list[i][j].replace("e-","D-")
    
    # Tabulate data
    cell_new.write(tabulate(line_list, tablefmt="plain", colalign=("left",)))
    cell_new.close()

#############################################################################
def modify_existing_cell(output_file, cell_file, ex, ey, ez):
    """This function adds displacements to cell files already 
    displaced by another field - increases error so not used.
    Out of date also (can be updated if needed).
    
    Parameters
    ----------
    output_file: str
        Path of CRYSTAL output file containing the unit cell data
    cell_file: str
        Path of initial cell file containing a form of unit cell
        data and charge data. Can be generated by CRYSTAL or user.
    ex: float
        Input for electric field in 'x' vector axis (Ex)
    ey: float
        Input for electric field in 'y' vector axis (Ey)
    ez: float
        Input for electric field in 'c' vector axis (Ez)
    
    Returns
    ----------
    None
        Writes to file, no returned output within python
    """
    
    displacements = convert_disp(output_file, ex, ey, ez)
    cell_old = open(cell_file).readlines()
    cell_new = []
    
    for i in range(len(cell_old)):
        cell_split = cell_old[i].split()
        for j in range(len(cell_split)):
            cell_split[j] = cell_split[j].replace("D+", "e+")\
                            .replace("D-", "e-")
        
        for j in range(3):
            cell_split[j+1] = str("{:.15e}".format(float(cell_split[j+1])
                              +float(displacements[3*i + j])))
        temp_list = []
        
        for j in range(len(cell_split)):
            temp_list.append(cell_split[j].replace("e+", "D+")\
            .replace("e-", "D-"))
        
        cell_new.append(temp_list)
    
    prev_field = ".".join(cell_file.split(".")[1:5])
    curr_field = "(%s,%s,%s)"%(str(ex), str(ey), str(ez))
    compound = cell_file.split(".")[0]
    new_cell_name = "%s.%s&%s.cell"%(compound, prev_field, curr_field)
    
    f = open(new_cell_name, 'w')
    f.write(tabulate(cell_new, tablefmt='plain', colalign=("left",)))
########################################################################

def cell_grid(born_file, hess_file, output_file, cell_file,
              ex_array, ey_array, ez_array, unit_source, *args, **kwargs):
    """This function generates a grid of displaced cell files
    in the target directory. Uses range of Ex, Ey, Ez values. As before,
    coordinates are in fractional units in the a,b,c basis.

    Parameters
    ----------
    born_file: str
        Path of born file
    hess_file: str
        Path of Hessian file
    output_file: str
        Path of CRYSTAL output file containing the unit cell data
    cell_file: str
        Path of initial cell file containing a form of unit cell
        data and charge data. Can be generated by CRYSTAL or user.
    ex_array: array
        Range of values for electric field input in 'x' vector axis (Ex)
    ey_array: array
        Range of values for electric field input in 'y' vector axis (Ey)
    ez_array: array
        Range of values for electric field input in 'z' vector axis (Ez)
    unit_source: str
        Three key words that determine the source of the initial
        non-disturbed unit cell fractional coordinates and the
        associated charge with each atom type
    *args: dict
        Charge data dictionary e.g {"Y":3.0,"MN":3.0,"O1":-2.0,
        "O2",-2.0}. This is the argument for when unit_source =
        'charge' is chosen, for manual user input
    
    Returns
    ----------
    directory: str
        Writes grid to target directory, output is the grid directory
        for crys2seward
    """

    # Generating new folder for results, named after ranges of E coordinates
    if len(ex_array) == 1 or len(ey_array) == 1 or len(ez_array) == 1:
        dirname = "Ex_%s_%s__Ey_%s_%s__Ez_%s_%s_cell"%\
                (str(round(ex_array[0], 3)),
                 str(round(ex_array[-1], 3)),str(round(ey_array[0], 3)),
                 str(round(ey_array[-1], 3)),str(round(ez_array[0], 3)),
                 str(round(ez_array[-1], 3)))

    else:
        dirname = "Ex_%s_%s_%s__Ey_%s_%s_%s__Ez_%s_%s_%s_cell"\
                %(str(round(ex_array[0], 3)),str(round(ex_array[-1], 3)),
                  str(round(abs(ex_array[0]-ex_array[1]), 3)),
                  str(round(ey_array[0], 3)), str(round(ey_array[-1], 3)),
                  str(round(abs(ey_array[0]-ey_array[1]), 3)),
                  str(round(ez_array[0], 3)), str(round(ez_array[-1], 3)),
                  str(round(abs(ez_array[0]-ez_array[1]), 3)))

    # Make new folder in target directory from modify_cell,
    # and print location
    out_dir = "/".join(output_file.split("/")[:-1])
    cell_dir = "/".join(cell_file.split("/")[:-1])

    # Determining where to write the files to
    # If inputs are paths, use them as the file destination in
    # priority order specified below
    if os.path.isdir(cell_dir):
        cd = cell_dir
        print("New cell files output to initial cell file directory")

    elif os.path.isdir(out_dir):
        cd = out_dir
        print("New cell files output to crystal output file directory")

    else:
        cd = os.getcwd()
        print("New cell files output to current working directory")

    directory = cd+"/"+dirname

    # Creating the target directory folder
    Path(directory).mkdir(parents=True, exist_ok=True)
 
    # Writing cell files
    for i in range(len(ex_array)):
        for j in range(len(ey_array)):
            for k in range(len(ez_array)):
                modify_cell(born_file, hess_file, output_file,
                            cell_file, dirname, ex_array[i], ey_array[j],
                            ez_array[k], unit_source, *args, **kwargs)

    print("Cell Grid written")
    return directory


def read_input(input_file):
    """Scans input file and calls cell_grid based on the contained info
    
    Parameters
    ----------
    input_file: str
        Path of user-generated input file containing required information
        to generate the cell file grid for a range of electric fields and
        displacements with cell_grid

    Returns
    ----------
    cd: str
        Cell grid written from input_file parsed info, output is the
        directory of this grid for crys2seward
    """
    
    # Parse for key words and take inputs for cell_grid
    inp = open(input_file).readlines()
    for i in range(len(inp)):
        inp_split = inp[i].split()
        
        if inp_split[0] == "born_file":
            born_file = "".join(inp_split[2:]).strip()
        
        elif inp_split[0] == "hess_file":
            hess_file = "".join(inp_split[2:]).strip()

        elif inp_split[0] == "crystal_file":
            crystal_file = "".join(inp_split[2:]).strip()
        
        elif inp_split[0] == "cell_init":
            cell_init = "".join(inp_split[2:]).strip()
        
        elif inp_split[0].lower() == "ex":
            ex_list = ast.literal_eval("".join([i for i in inp_split[2:]]))
        
        elif inp_split[0].lower() == "ey":
            ey_list = ast.literal_eval("".join([i for i in inp_split[2:]]))
        
        elif inp_split[0].lower() == "ez":
            ez_list = ast.literal_eval("".join([i for i in inp_split[2:]]))
        
        elif inp_split[0].lower() == "unit_source":
            unit_source = "".join(inp_split[2:]).strip()
        
        elif inp_split[0].lower() == "charge_dict":
            charge_dict = "".join(inp_split[2:]).strip()
            charge_dict = charge_dict.replace("'",'"')
            charge_dict = json.loads(charge_dict)
    
    # Zero array options
    if ex_list == [0, 0, 0]:
        ex_array = np.array([0])
    
    else:
        ex_array = np.arange(ex_list[0], ex_list[1]+ex_list[2], ex_list[2])
    
    if ey_list == [0, 0, 0]:
        ey_array = np.array([0])
    
    else:
        ey_array = np.arange(ey_list[0], ey_list[1]+ey_list[2], ey_list[2])

    if ez_list == [0, 0, 0]:
        ez_array = np.array([0])
    
    else:
        ez_array = np.arange(ez_list[0], ez_list[1]+ez_list[2], ez_list[2])

    # Charge unit_source requires separate option including charge_dict
    if unit_source == "charge":
        print("Cell grid is being generated...")
        cd = cell_grid(born_file, hess_file, crystal_file, cell_init,
                       ex_array, ey_array, ez_array, unit_source,
                       charge_dict)

    elif unit_source == "auto" or unit_source == "direct":
        print("Cell grid is being generated...")
        cd = cell_grid(born_file, hess_file, crystal_file, cell_init,
                       ex_array, ey_array, ez_array, unit_source)

    else:
        # Invalid unit_source again, cancel run
        sys.exit("unit_source is invalid, see documentation")

    return cd


def return_atoms(born_file, hess_file, output_file, cell_file, unit_source):
    """Lists atoms in the system to be used in the
    command line prompts

    Parameters
    ----------
    born_file: str
        Path of born file
    hess_file: str
        Path of Hessian file
    output_file: str
        Path to Crystal output file containing unit cell data
    cell_file: str
        Path to initial cell_file which can be used to generate coordinate
        and charge data for each atom
    unit_source: str
        Key word to determine source of unit cell information

    Returns
    ----------
    atoms: list
        Contains list of all atoms in the unit cell to be prompted for.
    """

    # Generate the atom dictionary to scan for atoms
    atom_dict = unit_cell(born_file, hess_file, output_file,
                          cell_file, 0, 0, 0, unit_source)[0]
    
    atoms = []
    
    # Append the atoms to the atoms list
    for group in atom_dict:
        for atom in atom_dict[group]:
            elem = atom_dict[group][atom]["element"]
            if elem not in atoms:
                atoms.append(elem)
    return atoms


def read_prompt():
    """Prompt entry of data required to call cell_grid

    No parameters or returns. Files are written and saved
    to a target directory
    """
    
    print("####################PROMPT INPUT#####################")
    print("Please enter file names if local to script directory,"+\
          "or full path if non-local")
    
    # Unit source prompt
    unit_source = input("Please give unit cell generator parameter"+\
                        " (direct,charge,auto): ").lower().split()[0]
    
    # If not recognised, end process
    if unit_source not in ["auto", "charge", "direct"]:
        sys.exit("Unit cell generator parameter not recognised")
    
    print("FILE PATHS: accepts crystal output or .DAT for Born and"+\
          " Hessian matrices")
    # Output and cell file path prompt
    born_file = input("Name/directory of born tensor source file: ")
    hess_file = input("Name/directory of Hessian matrix source file: ")
    crystal_file = input("Name/directory of crystal output file: ")
    cell_init = input("Name/directory of initial cell file containing"+\
                      " the unit cell data: ")
    
    # Electric field prompts in tuples
    print("Parameters for Ex array: ")
    ex_list = input("Please enter start,stop,step separated by commas"+\
                    " for Ex: ").split(",")
    
    print("Parameters for Ey array: ")
    ey_list = input("Please enter start,stop,step separated by commas"+\
                    " for Ey: ").split(",")
    
    print("Parameters for Ez array: ")
    ez_list = input("Please enter start,stop,step separated by commas"+\
                    " for Ez: ").split(",")
    
    # If zero list, make zero for whole grid
    if ex_list == ["0", "0", "0"]:
        ex_array = np.array([0])
    
    # Else generate input arrays
    else:
        ex_array = np.arange(float(ex_list[0]),
                             float(ex_list[1])+float(ex_list[2]),
                             float(ex_list[2]))
    
    if ey_list == ["0", "0", "0"]:
        ey_array = np.array([0])
    
    else:
        ey_array = np.arange(float(ey_list[0]),
                             float(ey_list[1])+float(ey_list[2]),
                             float(ey_list[2]))
    
    if ez_list == ["0", "0", "0"]:
        ez_array = np.array([0])
    
    else:
        ez_array = np.arange(float(ez_list[0]),
                             float(ez_list[1])+float(ez_list[2]),
                             float(ez_list[2]))
    
    # Asking for charge data for specific atoms
    if unit_source == "charge":
        element_charge = {}
        atoms = return_atoms(crystal_file, cell_init, unit_source)
        print("charge parameter has been chosen, please center charge data")
        for i in range(len(atoms)):
            element_charge[atoms[i]] = float(input("Please enter the float"+\
                    " value of the charge for %s: "%atoms[i]))
        
        print("###########################END OF INPUT###################"+\
                "##########")
        print("Cell grid is being generated...")
        
        # Generate grid
        cell_grid(born_file, hess_file, crystal_file, cell_init,
                  ex_array, ey_array, ez_array, unit_source,
                  element_charge)
    
    else:
        print("###########################END OF INPUT##################"+\
                "###########")
        print("Cell grid is being generated...")
        
        # Generate grid
        cell_grid(born_file, hess_file, crystal_file, cell_init,
                  ex_array, ey_array, ez_array, unit_source)


def input_or_prompt(crys2sew_bool):
    """Logic for receiving input file or prompting user
    
    Parameters
    ----------
    crys2sew_bool: bool
        True: using python script as package, don't want to run
        anything by default, stop the function
        False: Python script being used directly, run either the
        prompt or parsing of the third input file

    Returns
    ----------
    None
    """

    # If input bool is True, functions are called by crys2seward.py
    if crys2sew_bool:
        pass

    # If an input file argument is passed on the command line
    # run read_input
    elif len(sys.argv) >= 2:
        read_input(sys.argv[1])
    
    # Else prompt the user with read_prompt
    else:
        read_prompt()

#######################################################################
def basis_set_conv(born_file, hess_file, output, output2):
    """Generates the conversion matrix between two systems of the same
    coordinate basis but different basis sets and applies it to the
    corresponding Hessian.

    Parameters
    ----------
    born_file: str
        Path of born file
    hess_file: str
        Path of Hessian file
    output_file: str
        Path to Crystal output file containing unit cell data

    Returns
    ----------
    tensor_hess
    """

    prim2cart_matrix = np.zeros((3, 3))
    prim2cart_matrix2 = np.zeros((3, 3))
    hess_input(born_file, hess_file, output)

    with open(output, 'r') as file:
        crystal = file.readlines()
 
    with open(output2, 'r') as file:
        crystal2 = file.readlines()

    for i in range(len(crystal)):
        if crystal[i].split() == ['DIRECT', 'LATTICE', 'VECTORS',
                                  'CARTESIAN', 'COMPONENTS', '(ANGSTROM)']:
            start = i + 2
            break

    conv_basis_data = crystal[start : start+3]

    for i in range(len(crystal2)):
        if crystal2[i].split() == ['DIRECT', 'LATTICE', 'VECTORS',
                                  'CARTESIAN', 'COMPONENTS', '(ANGSTROM)']:
            start2 = i + 2
            break

    conv_basis_data2 = crystal2[start2 : start2+3]

    for i in range(len(conv_basis_data)):
        cbd_split = conv_basis_data[i].split()
        for j in range(len(cbd_split)):
            prim2cart_matrix[i, j] = float(cbd_split[j])

    for i in range(len(conv_basis_data2)):
        cbd_split2 = conv_basis_data2[i].split()
        for j in range(len(cbd_split2)):
            prim2cart_matrix2[i, j] = float(cbd_split2[j])

    prim2cat_matrix3 = np.matmul(np.linalg.inv(prim2cart_matrix),
                                 prim2cart_matrix2)
######################################################################

# Test inputs for functions and calling input_or_prompt

# input_file = "/home/joshhorswill10/Documents/git_new/joshua_3/Examples/"+\
#              "crys2seward_examples/ht.frequence.B1PW_PtBs.loto.out"
# cell_first = "/home/joshhorswill10/Documents/git_new/joshua_3/Examples/"+\
#               crys2seward_examples/ymno3.new.cell"
# ex, ey, ez = np.arange(0,0.6,.1), np.arange(0,0.21,0.01),\
#              np.arange(0,1.2,0.2)
# {'Y':3.0, 'MN':3.0, 'O1':-2.0, 'O2':-2.0}
# b_dat = "YMnO3/BORN_B1Pw_loto.DAT"
# h_dat = "YMnO3/HESSIEN.DAT"
# crys_out = "../frequence.B1PW.loto.out"
# CuO_P1 = "/home/joshhorswill10/Documents/git_new/joshua_3/Examples/"+\
#          "crys2seward_examples/CuO_examples/CuO-P1.AFM.Phonons.156273.out"
# CuO_C2 = "/home/joshhorswill10/Documents/git_new/joshua_3/Examples/"+\
#          "crys2seward_examples/CuO_examples/CuO-C2.AFM.Phonons.155494.out"
# cell_CuO = "/home/joshhorswill10/Documents/git_new/joshua_3/Examples/"+\
#            "crys2seward_examples/CuO_examples/CuO_P1.cell"
# print("Cartesian displacement matrices: ")
# solve_equation(CuO_P1, CuO_P1, CuO_P1, 1e7, 0, 0)

# print("Converted displacements only: ")
# print(convert_disp(CuO_P1, CuO_P1, CuO_P1, 1e7, 0, 0)[0])
# basis_set_conv(CuO_P1, CuO_P1, CuO_P1, CuO_C2)
# print(unit_cell(CuO_P1, CuO_P1, CuO_P1, cell_CuO, 0, 0, 0, "direct")[0])

input_or_prompt(True)
