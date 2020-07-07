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
    atom_array = np.arange(1,atom_num+1,1)
    matrix_dict = {"atom%d"%i:np.zeros((3,3)) for i in atom_array}
    born = open(born_file).readlines()
    for i in range(len(born)):
        born_split = born[i].split()
        if born_split[:3] == ['BORN', 'CHARGE', 'TENSOR.']:
            start = i+1
            end = i+1 + (7*atom_num)
            break
    tensor_data = born[start:end]
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
        for j in range(len(hess_split)):
            if float(hess_split[j])<= 1e-15:
                hess_data.append(0.0)
            else:
                hess_data.append((float(hess_split[j]),i*len(hess_split)+j))
    hess_copy = copy.deepcopy(hess_data) 
    hess_data = [i for i in hess_data if i != 0.0]
    hess_pairs = []
    hess_diag = []
    flag = True
    while flag == True:
        #print("initial element length ",len(hess_data),"twice pair length + diag length ",2*len(hess_pairs)+len(hess_diag))
        if len(hess_data) == 1:
            hess_diag.append(hess_data[0])
            hess_data.remove(hess_data[0])
            break
        if len(hess_data) == 0:
            flag = False
            break
        for i in range(1,len(hess_data)):
            if abs(hess_data[0][0]-hess_data[i][0]) <= thresh:
                mean = (hess_data[0][0]+hess_data[i][0])/2.0
                hess_pairs.append([mean,hess_data[0][1],hess_data[i][1]])
                rem1 = hess_data[0]
                rem2 = hess_data[i]
                hess_data.remove(rem1),hess_data.remove(rem2)
                break
            elif i == len(hess_data)-1:
                hess_diag.append([hess_data[0][0],hess_data[0][1]])
                hess_data.remove(hess_data[0])
    return hess_pairs,hess_diag,hess_copy

def hessian(born_file,hess_file,thresh):
    "Returns Hessian data as a square matrix (unfinished)"
    m_dim = 3*num_atoms(born_file)
    tensor_hess = np.zeros((m_dim,m_dim))
    hess_info = parse_hessian(hess_file,thresh)
    hess = open(hess_file).readlines()
    #Waiting on example

def solve_equation(born_file,hess_file,thresh,ex,ey,ez):
    "Using numpy.linalg.solve for an ax = b equation"
    #q_tensor = born_tensor(born_file)
    #h_tensor = hessian(born_file,hess_file,thresh)
    #e_vector = np.ones((len(h_tensor)))
    #Placeholder matrix so we can keep writing the script
    m_dim = 192
    h_tensor = np.random.rand(m_dim,m_dim)
    q_tensor = np.random.rand(m_dim,m_dim)
    e_vector = np.ones((m_dim))
    for i in range(len(e_vector)):
        if i % 3 == 0:
            e_vector[i] = ex
        elif i% 3 == 1:
            e_vector[i] = ey
        else:
            e_vector[i] = ez
    a = -1*h_tensor
    b = q_tensor.dot(e_vector)
    displacement = np.linalg.solve(a,b)
    print("Negative Hessian: ")
    print(a)
    print("Born tensor: ")
    print(q_tensor)
    print("Electric field: ")
    print(e_vector)
    print("Displacement solutions: ")
    print(displacement)

def print_indexes(hess_file,thresh):
    "Prints all paired and unique values with their frequencies"
    print("Pair indexes")
    print([i[1:] for i in parse_hessian(hess_file,thresh)[0]])
    print("Diag indexes")
    print([i[1] for i in parse_hessian(hess_file,thresh)[1]])

def print_diag(hess_file,thresh,rev):
    """Prints diagonal elements in either ascending
    or descending order depending on rev boolean"""
    with open("diagonal.txt",'w') as diag:
        sort_diag = parse_hessian(hess_file,thresh)[1]
        sort_diag = [i[0] for i in sort_diag]
        sort_diag.sort(reverse=rev)
        print("Length ", len(sort_diag))
        for item in sort_diag:
            diag.write("%s\n"%str(item))

def num_diag(hess_file):
    thresh_range = np.ones((15))
    for i in range(len(thresh_range)):
        if i == 0:
            thresh_range[i] = 1e-3
        else:
            thresh_range[i] = 1e-1*thresh_range[i-1]
    diag_dict = {str(i):None for i in thresh_range}
    for i in range(len(thresh_range)):
        diag_dict['%s'%str(thresh_range[i])] = len(parse_hessian("Pm.rePhonons.2123667.HESSFREQ.DAT",thresh_range[i])[1])
    print("Format is 'threshold':number of diagonal elements")
    print(diag_dict)

num_diag("Pm.rePhonons.2123667.HESSFREQ.DAT")
#solve_equation('frequence.B1PW_PtBs.loto.out',"Pm.rePhonons.2123667.HESSFREQ.DAT",1,2,3)
#print_diag("Pm.rePhonons.2123667.HESSFREQ.DAT",1e-4,True)
#print_indexes("Pm.rePhonons.2123667.HESSFREQ.DAT")
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
