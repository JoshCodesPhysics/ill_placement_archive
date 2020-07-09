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
            tensor[i,j+start_column] = m_dict['atom%d'%atom_num][i%3,j]
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

def hessian(born_file,hess_file):
    "Returns L Hessian matrix as a square matrix"
    #Determines dimensions of Hessian
    #m_dim = 3*num_atoms(born_file)
    m_dim = 192
    #Empty square matrix
    tensor_hess = np.zeros((m_dim,m_dim))
    hess = open(hess_file).readlines()
    final_block = False
    #Finding relevant data range
    for i in range(len(hess)):
        hess_split = hess[i].split()
        if hess_split == [str(i) for i in range(1,len(hess_split)+1)]\
                and len(hess_split) > 0 and hess[i-1] == '\n':
            start = i
            print("start = ",i)
        if ((str(m_dim) in hess_split and str(m_dim - 1) in hess_split)\
                or (str(m_dim) in hess_split and str(m_dim+1) in hess_split)\
                or (hess_split == [str(m_dim)])):
            final_block = True
            print("Final block = ",final_block)
        if len(hess_split) > 0:
            if hess_split[0] == str(m_dim+3) and final_block == True:
                end = i+1
                print("end = ", end)
    #defining new range
    hess_data = hess[start:end]
    print(hess_data[0],hess_data[-1])
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
            if i <= j:
                pass
            else:
                tensor_hess[j,i] = tensor_hess[i,j]
    #Print example of symmetries
    """entries = 10
    for i in range(entries):
        for j in range(entries):
            if i != j:
                print(tensor_hess[i,j],tensor_hess[j,i],(i,j),(j,i))
    """
    #Add [0] to function for L matrix and [1] for symmetric matrix
    return l_matrix,tensor_hess

def solve_equation(born_file,hess_file,ex,ey,ez):
    "Using numpy.linalg.solve for an ax = b equation"
    #q_tensor = born_tensor(born_file)
    h_tensor = hessian(born_file,hess_file)[1]
    m_dim = len(h_tensor)
    e_vector = np.ones((m_dim))
    #Placeholder matrix so we can keep writing the script   
    q_tensor = np.random.rand(m_dim,m_dim)
    #Generating e_vector from inputs
    for i in range(len(e_vector)):
        if i % 3 == 0:
            e_vector[i] = ex
        elif i% 3 == 1:
            e_vector[i] = ey
        else:
            e_vector[i] = ez
    #Generating matrices for linalg solve
    a = -1*h_tensor
    b = q_tensor.dot(e_vector)
    displacement = np.linalg.solve(a,b)
    #Printing components of linear equation including the solution
    print("Negative Hessian: ")
    print(a)
    print("Born tensor: ")
    print(q_tensor)
    print("Electric field: ")
    print(e_vector)
    print("Displacement solutions: ")
    print(displacement)

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
        diag_dict['%s'%str(thresh_range[i])] = len(parse_hessian("Pm.rePhonons.2123667.HESSFREQ.DAT",thresh_range[i])[1])
    print("Format is 'threshold':number of diagonal elements")
    print(diag_dict)

#solve_equation('frequence.B1PW_PtBs.loto.out','Pm.reOptGeom.2164102.out',2,0,0)
print(hessian('frequence.B1PW_PtBs.loto.out','Pm.rePhonons.7094721.out')[1])
#num_diag("Pm.rePhonons.2123667.HESSFREQ.DAT",5,1e-3)
#print_diag("Pm.rePhonons.2123667.HESSFREQ.DAT",1e-12,True)
#print_indexes("Pm.rePhonons.2123667.HESSFREQ.DAT",1e-4)
#hessian("Pcab.reOptGeom.364624.out","Pcab.reOptGeom.364624.HESSOPT.DAT")
#tensor = born_tensor('frequence.B1PW_PtBs.loto.out')
#print("Born tensor: ",tensor)
#print("#######################")
#print("tensor*inverse: ",tensor.dot(np.linalg.inv(tensor)))

"""def cmd_input_or_prompt():
    if len(sys.argv) >= 2:
       pass
    else:
        pass
"""
