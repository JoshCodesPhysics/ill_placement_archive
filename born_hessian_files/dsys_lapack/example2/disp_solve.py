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
np.set_printoptions(precision=15)

def num_atoms(born_file):
    "Determines number of atoms in the unit cell from parse"
    born = open(born_file).readlines()
    for i in range(len(born)):
        born_split = born[i].split()
        if born_split[:5] == ['ATOMS', 'IN', 'THE', 'ASYMMETRIC', 'UNIT']:
            atom_num = int(born_split[-1])
            break
    return atom_num

def matrix_data(born_file):
    "Stores matrix data for each atom as 3x3 matrix corresponding to the atom key"
    atom_num = num_atoms(born_file)
    #Generating atom data and empty storage matrix
    atom_array = np.arange(1,atom_num+1,1)
    matrix_dict = {"atom%d"%i:np.zeros((3,3)) for i in atom_array}
    born = open(born_file).readlines()
    #Searching for title and truncating relevant data
    for i in range(len(born)):
        born_split = born[i].split()
        if born_split[:3] == ['BORN', 'CHARGE', 'TENSOR.']:
            start = i+1
            end = i+1 + (7*atom_num)
            break
    tensor_data = born[start:end]
    #Appending data to dictionary as 3x3 matrices
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
                        if abs(float(tensor_data[m_index].split()[k+1])) < 1e-12:
                            matrix_dict["atom%d"%atom_index][j,k] = 0e0
                        else:
                            matrix_dict["atom%d"%atom_index][j,k] = float(tensor_data[m_index].split()[k+1])
    return matrix_dict

def born_tensor(born_file):
    "Appends this matrix data to a larger matrix where each 3x3 matrix is a diagonal element of the larger matrix"
    m_dict = matrix_data(born_file)
    m_dim = 3*num_atoms(born_file)
    tensor = np.zeros((m_dim,m_dim))
    start_column = 0
    atom_num = 1
    #Looping through tensor columns and adding stepped matrices along the diagonal
    for i in range(m_dim):
        if i != 0 and i%3 ==0:
            start_column += 3
            atom_num += 1
        for j in range(3):
            tensor[i,j+start_column] = float(m_dict['atom%d'%atom_num][i%3,j])
    return tensor

def parse_hessian(hess_file,thresh):
    "Parses Hessian to find index of pairs and diagonal elements"
    hess_data = []
    hess = open(hess_file).readlines()
    for i in range(len(hess)):
        hess_split = hess[i].split()
        #Appending all values to hess_data from input file, unless the value is very small, then it is appended as zero
        for j in range(len(hess_split)):
            if abs(float(hess_split[j]))<= 1e-15:
                hess_data.append(0.0)
            else:
                #Appending data and keeping the original index for later
                hess_data.append((float(hess_split[j]),i*len(hess_split)+j))
    #Keeping a copy of the data with zeros included
    hess_copy = copy.deepcopy(hess_data)
    #Removing zeros
    hess_data = [i for i in hess_data if i != 0.0]
    hess_pairs = []
    hess_diag = []
    flag = True
    #While there are still elements in hess_data, append and remove pairs
    #with their respective indexes as well as unique elements with their indexes.
    while flag == True:
        #Uncomment below statement if you want a 'progress bar'
        #print("initial element length ",len(hess_data),"twice pair length + diag length ",2*len(hess_pairs)+len(hess_diag))
        if len(hess_data) == 1:
            #Removes bug for last element
            hess_diag.append(hess_data[0])
            hess_data.remove(hess_data[0])
            break
        #Break once list is empty
        if len(hess_data) == 0:
            flag = False
            break
        for i in range(1,len(hess_data)):
            #If the difference between two elements satisfies the threshold
            #Add them and their indexes to a pair list and remove both from hess_data
            if abs(hess_data[0][0]-hess_data[i][0]) <= thresh:
                mean = (hess_data[0][0]+hess_data[i][0])/2.0
                hess_pairs.append([mean,hess_data[0][1],hess_data[i][1]])
                rem1 = hess_data[0]
                rem2 = hess_data[i]
                hess_data.remove(rem1),hess_data.remove(rem2)
                break
            #If there are no pairs for this value, add it to diagonal list with it's index and remove
            elif i == len(hess_data)-1:
                hess_diag.append([hess_data[0][0],hess_data[0][1]])
                hess_data.remove(hess_data[0])
    #return triple of pairs, diagonal and original data lists
    return hess_pairs,hess_diag,hess_copy

def print_indexes(hess_file,thresh):
    "Prints all paired and unique values with their indexes"
    print("Pair indexes")
    print([i[1:] for i in parse_hessian(hess_file,thresh)[0]])
    print("Diag indexes")
    print([i[1] for i in parse_hessian(hess_file,thresh)[1]])

def print_diag(hess_file,thresh,rev):
    """Writes diagonal elements to a text file in either ascending
    or descending order depending on rev boolean"""
    with open("diagonal.txt",'w') as diag:
        sort_diag = parse_hessian(hess_file,thresh)[1]
        sort_diag = [i[0] for i in sort_diag]
        #Sorting ascending or descending, rev = True creates descending
        sort_diag.sort(reverse=rev)
        print("Length ", len(sort_diag))
        for item in sort_diag:
            diag.write("%s\n"%str(item))

def num_diag(hess_file,nthresh,threshstart):
    "Prints a dictionary containing number of diagonal elements for a range of thresholds"
    "nthresh gives number of thresholds to test, and threshstart gives starting value"
    thresh_range = np.ones((nthresh))
    #Generating threshold values
    for i in range(len(thresh_range)):
        if i == 0:
            thresh_range[i] = threshstart
        else:
            #One order of magnitude between each index
            thresh_range[i] = 1e-1*thresh_range[i-1]
    #Empty dictionary for writing
    diag_dict = {str(i):None for i in thresh_range}
    #Loop writing all diagonal list lengths for given thresholds
    for i in range(len(thresh_range)):
        diag_dict['%s'%str(thresh_range[i])] =\
                len(parse_hessian("Pm.rePhonons.2123667.HESSFREQ.DAT",thresh_range[i])[1])
    print("Format is 'threshold':number of diagonal elements")
    print(diag_dict)


def hessian(born_file,hess_file):
    "Returns L Hessian matrix as a square matrix"
    #Determines dimensions of Hessian
    m_dim = 3*num_atoms(born_file)
    #Empty square matrix
    tensor_hess = np.zeros((m_dim,m_dim))
    hess = open(hess_file).readlines()
    final_block = False
    start_found = False
    #Finding relevant data range
    for i in range(len(hess)):
        hess_split = hess[i].split()
        #Start condition for Pm.rePhonons file
        if hess_split == ['LOW', 'HALF', 'OF', 'HRED', 'AFTER', 'SYMMETRISATION']:
            start = i+4
            start_found = True
        #Start condition for Pm.reOptGeom file
        elif hess_split == [str(i) for i in range(1,len(hess_split)+1)]\
                and len(hess_split) > 0 and hess[i-1] == '\n'\
                and len(hess[i+1].split()) == 2 and start_found == False:
            start = i
            start_found = True
        #Detecting final block to prepare for end truncation
        if ((str(m_dim) in hess_split and str(m_dim - 1) in hess_split)\
                or (str(m_dim) in hess_split and str(m_dim+1) in hess_split)\
                or (hess_split == [str(m_dim)])):
            if start_found == True:
                final_block = True
        #End detected, break the loop
        if len(hess_split) > 0:
            if hess_split[0] == str(m_dim) and final_block == True and start_found == True:
                end = i+1
                break
    #defining new range
    hess_data = hess[start:end]
    #Generalised column numbers for multiple format possibilities
    ncolumns = len(hess_data[0].split())
    #First column starts at 0, increase by number of columns for every data block
    first_column = 0
    #If we are in the last block and we hit the last three rows and columns - break
    last_block = False
    for i in range(len(hess_data)):
        #Splitting lines into individual elements with no spacing
        rowsplit = hess_data[i].split()
        #print("First column ",first_column)
        #If we hit an empty line, increase the column number by ncolumns
        if len(rowsplit) == 0:
            first_column += ncolumns
            #print("empty space")
        #If we hit a row of column labels, do nothing except for if it is the last block
        elif rowsplit[:2] == [str(first_column+1),str(first_column+2)]:
            #print("Column label row")
            if str(int(m_dim+1)) in rowsplit:
                last_block = True
        #Data lines
        else:
            #Row number is first element given
            row = int(rowsplit[0])-1
            #If we hit last block and last three rows, stop the loop
            if last_block == True and row >= m_dim:
                break
            #If we hit the last three rows in a non-final block, just don't take the data
            elif last_block == False and row >= m_dim:
                pass
            #Appending values to tensor
            else:
                for j in range(1,len(rowsplit)):
                    tensor_hess[row,j+first_column-1] = float(rowsplit[j])
            #print("Appending data for rows %d-%d"%(first_column,first_column+ncolumns),"Row: ",row)
    #Copy of L matrix if it is needed
    l_matrix = copy.deepcopy(tensor_hess)
    #Symmetrising
    for i in range(m_dim):
        for j in range(m_dim):
            if i >= j:
                if abs(tensor_hess[i,j]) <= 1e-12:
                    tensor_hess[i,j] = 0
            else:
                if abs(tensor_hess[j,i]) <= 1e-12:
                    tensor_hess[i,j] = 0
                else:
                    tensor_hess[i,j] = float(tensor_hess[j,i])
    #Print example of symmetries
    """entries = 10
    for i in range(entries):
        for j in range(entries):
            if i != j:
                print(tensor_hess[i,j],tensor_hess[j,i],(i,j),(j,i))
    """
    #Add [0] suffix to function for L matrix and [1] for symmetric matrix
    return l_matrix,tensor_hess

def convert_coords(output_file):
    """Finds conversion matrix and transposes it. Also parses for
    fractional coordinate data"""
    #Allocating number of unit cell atoms and conversion matrices
    atom_num = num_atoms(output_file)
    out = open(output_file).readlines()
    conv_matrix = np.zeros((3,3))
    conv_inverse= np.zeros((3,3))
    #Angstrom to Bohr conversion factor
    ang2at = float(0.529177210903)
    #Parsing for data strips
    for i in range(len(out)):
        out_split = out[i].split()
        if out_split == ['DIRECT', 'LATTICE', 'VECTORS', 'CARTESIAN',\
                'COMPONENTS', '(ANGSTROM)']:
            start = i+2
            end = start+3
        elif out_split == ['COORDINATES', 'OF', 'THE', 'EQUIVALENT',\
                'ATOMS', '(FRACTIONARY', 'UNITS)']:
            start2 =i+4
        elif out_split[:5] == ['NUMBER','OF','SYMMETRY','OPERATORS',':']:
            end2 = i-1
        elif out_split == ['*', 'ATOM', 'X(ANGSTROM)', 'Y(ANGSTROM)',\
                'Z(ANGSTROM)']:
            start3 = i+2
            end3 = start3+atom_num
            break
    conv_data = out[start:end]
    frac_data = out[start2:end2]
    test_data = out[start3:end3]
    #Appending conversion matrix
    for i in range(len(conv_data)):
        conv_split = conv_data[i].split()
        for j in range(len(conv_split)):
            conv_matrix[i,j] = ang2at*float(conv_split[j])
    #Transposing it so we can solve for u,v,w
    conv_matrix = np.transpose(conv_matrix)
    conv_inverse = np.linalg.inv(conv_matrix)
    """
    #print("Translation matrix before being converted from Angstrom to Bohr: ")
    #print(conv_matrix/ang2at)
    frac_coords = np.zeros((atom_num,3))
    row = 0
    for i in range(len(frac_data)):
        frac_split = frac_data[i].split()
        if len(frac_split) > 0:
            for j in range(5,len(frac_split)):
                frac_coords[row,j-5] = frac_split[j]
            row += 1
    #print("Desired coordinates:")
    #print(frac_coords)
    
    test_coords = np.zeros((atom_num,3))
    for i in range(len(test_data)):
        test_split = test_data[i].split()
        for j in range(3,len(test_split)):
            if abs(float(test_split[j])) <= 1e-12:
                test_coords[i,j-3] = 0.0
            else:
                test_coords[i,j-3] = ang2at*float(test_split[j])
    print("########################################")
    print("Cartesian primitive cell coordinates before Bohr conversion:")
    print(test_coords/ang2at)
    print("########################################") 
    test_coords2 = np.zeros((atom_num,3))
    test_coords3 = np.zeros((atom_num,3))
    for i in range(len(test_coords)):
        xyz = test_coords[i]
        conv_coords = np.linalg.solve(conv_matrix,xyz)
        #conv_coords = np.matmul(conv_inverse,xyz)
        for j in range(len(conv_coords)):
            if abs(conv_coords[j]) > 0.5:
                conv_coords[j] = conv_coords[j] % 1
        test_coords2[i] = conv_coords
    for i in range(len(frac_coords)):
        abc = frac_coords[i]
        orig_coords = np.matmul(np.transpose(conv_matrix),abc)
        test_coords3[i] = orig_coords
    print("Final coordinates: ")
    print(test_coords2)
    print("#######################################")
    print("Original reconverted coordinates: ")
    print(test_coords3/ang2at)
    """
    return [conv_matrix,frac_data]

def convert_coords2(out_file):
    """Parses for lattice parameters and converts them into unit cell edge
    vectors in cartesian coordinates"""
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
    a_vector = np.array([a,0.0,0.0])
    b_vector = np.array([b*np.cos(gamma),b*np.sin(gamma),0])
    omega = a*b*c*np.sqrt(1-np.cos(alpha)*np.cos(alpha)-np.cos(beta)*np.cos(beta)\
            -np.cos(gamma)*np.cos(gamma)\
            +2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
    c_vector = np.array([c*np.cos(beta),\
            c*((np.cos(alpha)-np.cos(beta)*np.cos(gamma))/(np.sin(gamma))),\
            omega/(a*b*np.sin(gamma))])
    lattice_vectors = np.array([a_vector,b_vector,c_vector])
    unit_vectors = np.array([float(vec/np.sqrt(vec.dot(vec))) for vec in lattice_vectors])
    return lattice_vectors, unit_vectors

def save_hess_born(born_file,hess_file,ex,ey,ez):
    "Writes arrays to binary files for fortran processing"
    "Enter electric field components in kV/m"
    h_tensor = hessian(born_file,hess_file)[1]
    q_tensor = born_tensor(born_file)
    e_vector = np.zeros((len(h_tensor)))
    conversion = float(1.944690381e-9)
    #Generating e_vector from inputs
    for i in range(len(e_vector)):
        if i % 3 == 0:
            e_vector[i] = float(conversion*ex)
        elif i% 3 == 1:
            e_vector[i] = float(conversion*ey)
        else:
            e_vector[i] = float(conversion*ez)
    dim_array = np.array([len(h_tensor)])

    h_tensor.tofile('hess_matrix')
    q_tensor.tofile('born_matrix')
    e_vector.tofile('e_vector')
    dim_array.tofile('n_array')
    #np.save('hess_matrix',h_tensor,allow_pickle=True,fix_imports=False)
    #np.save('born_matrix',q_tensor,allow_pickle=True,fix_imports=False)
    print("q,E,H matrices written to binary files")

def solve_equation(born_file,hess_file,ex,ey,ez):
    "Using numpy.linalg.solve for an ax = b equation"
    "Give electric field in terms of kV/m and it will"
    "Be converted into atomic electric field units"
    q_tensor = born_tensor(born_file)
    h_tensor = hessian(born_file,hess_file)[1]
    m_dim = len(h_tensor)
    e_vector = np.zeros((m_dim))
    conversion = float(1.94469038e-9)
    #Placeholder matrix so we can keep writing the script   
    #q_tensor = np.random.rand(m_dim,m_dim)
    #Generating e_vector from inputs
    for i in range(len(e_vector)):
        if i % 3 == 0:
            e_vector[i] = float(conversion*ex)
        elif i% 3 == 1:
            e_vector[i] = float(conversion*ey)
        else:
            e_vector[i] = float(conversion*ez)
    #Generating matrices for linalg solve
    a = float(-1.0)*h_tensor
    b = np.matmul(q_tensor,e_vector)
    displacement = np.linalg.solve(a,b)
    #save_hess_born(born_file,hess_file,ex,ey,ez)
    #Printing components of linear equation including the solution
    #print("Negative Hessian: ")
    #print(a)
    #print("B matrix: ")
    #print(b)
    #print("Born tensor: ")
    #print(q_tensor)
    #print("Electric field: ")
    #print(e_vector)
    #print("Displacement solutions: ")
    #print(displacement)
    #print("Error matrix: ")
    #print(abs(np.matmul(a,displacement)-b))
    return displacement

def convert_disp(output_file,ex,ey,ez):
    displacement = solve_equation(output_file,output_file,ex,ey,ez)
    conv_data = convert_coords(output_file)
    conv_matrix = conv_data[0]
    frac = conv_data[1]
    convdisp = []
    for i in range(len(displacement)):
        if i == len(displacement)-1:
            pass
        elif i % 3 == 0:
            xyz = np.array(displacement[i:i+3])
            newcoords = np.linalg.solve(conv_matrix,xyz)
            for i in range(len(newcoords)):
                convdisp.append(newcoords[i])
    convdisp = np.array(convdisp)
    #print("Converted displacements")
    #print(convdisp)
    return convdisp,frac

def unit_cell(output_file,cell_file,ex,ey,ez):
    cell = open(cell_file).readlines()
    atom_num = num_atoms(output_file)
    conv_data = convert_disp(output_file,ex,ey,ez)
    convdisp = conv_data[0]
    frac_data = conv_data[1]
    num_groups = len(frac_data)-atom_num+1
    group_nums = [i for i in range(1,num_groups+1)]
    atom_nums = [i for i in range(1,atom_num+1)]
    atom_list = []
    group_list = []
    for i in atom_nums:
        atom_list.append("atom %d"%i)
    for i in group_nums:
        group_list.append("group %d"%i)
    gap_index = []
    for i in range(len(frac_data)):
        frac_split = frac_data[i].split()
        if len(frac_split) == 0:
            gap_index.append(i)
    #print(gap_index)
    groups = []
    for i in range(len(gap_index)):
        if i == 0:
            groups.append(frac_data[:gap_index[i]])
        elif i == len(gap_index)-1:
            groups.append(frac_data[gap_index[i-1]+1:gap_index[i]])
            groups.append(frac_data[gap_index[i]+1:])
        else:
            groups.append(frac_data[gap_index[i-1]+1:gap_index[i]])
    element_charge = {}
    for i in range(len(cell)):
        cell_split = cell[i].split()
        if cell_split[0][1].isdigit():
            element = cell_split[0][0].upper()
        else:
            element = cell_split[0].upper()
        charge = float(cell_split[-1])
        if (element,charge) not in element_charge.items():
            element_charge[element]=charge
    atom_dict = {key:{} for key in group_list}
    parsed_atoms = {}
    dist_atoms = {}
    for i in range(len(groups)):
        for j in range(len(groups[i])):
            coords_split = groups[i][j].split()
            #print(i,j,coords_split)
            atom = coords_split[4]
            atom_dict["group %d"%(i+1)]["atom %d"%(int(coords_split[0]))]=\
                        {"coords":[float(coords_split[5+k]) for k in range(3)],\
                        "element":atom,"charge":element_charge[atom]}
    for i in range(len(groups)):
        for j in range(len(groups[i])):
            coords_split = groups[i][j].split()
            atom = coords_split[4]
            if atom not in parsed_atoms and j == 0:
                    parsed_atoms[atom] = i+1
            elif coords_split[4] in parsed_atoms:
                if j > 0:
                    pass
                elif atom not in dist_atoms:
                    dist_atoms[atom] = 2
                    #print("Group of atom before: ",parsed_atoms[atom])
                    prev_dict = atom_dict["group %d"%parsed_atoms[atom]]
                    for k in prev_dict:
                        #print("prev_dict being written for %s"%k)
                        prev_dict[k]["element"] = "%s1"%atom
                    current_dict = atom_dict["group %d"%(i+1)]
                    #print("group of current atom: ",i+1)
                    for k in current_dict:
                        #print("current_dict being written for %s"%k)
                        current_dict[k]["element"] = "%s2"%atom
                else:
                    dist_atoms[atom] += 1
                    current_dict = atom_dict["group %d"%(i+1)]
                    for k in current_dict:
                        current_dict[k]["element"] = "%s%d"%(atom,dist_atoms[atoms])
    #print(atom_dict)
    return atom_dict,convdisp

def modify_cell(output_file,cell_file,ex,ey,ez):
    unit_data = unit_cell(output_file,cell_file,ex,ey,ez)
    atom_dict = unit_data[0]
    convdisp = unit_data[1]
    fsplit = cell_file.split('.')
    fsplit[1:1] = ['new']
    fnew = '.'.join(fsplit)
    cell_new = open(fnew,'w+')
    line_list = []
    for group in atom_dict:
        for atom in atom_dict[group]:
            temp_dict = atom_dict[group][atom]
            writelist = [str(temp_dict["element"]),None,None,None,str("{:.12e}".format(temp_dict["charge"]))+'\n']
            for i in range(3):
                writelist[i+1] = str("{:.12e}".format(temp_dict["coords"][i]))
            writeline = "     ".join(writelist)
            line_list.append(writeline.split())

    for i in range(len(line_list)):
        for j in range(1,len(line_list[i])):
            line_list[i][j] = line_list[i][j].replace("e+","D+")
            line_list[i][j] = line_list[i][j].replace("e-","D-")
    cell_new.write(tabulate(line_list,tablefmt="plain",colalign=("left",)))

input_file = "ht.frequence.B1PW_PtBs.loto.out"
cell_file = "ymno3.cell"
#solve_equation(input_file,input_file,.1,.1,.1)
#save_hess_born(input_file,input_file,.1,.1,.1)
#convert_coords(input_file)
#convert_disp(input_file,.1,.1,.1)
modify_cell(input_file,cell_file,.1,.1,.1)
