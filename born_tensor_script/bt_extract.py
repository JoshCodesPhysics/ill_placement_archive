#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 08:35:02 2020

@author: joshhorswill10
"""
import numpy as np
import sys
import itertools as it

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

def hessian(born_file,hess_file):
    "Returns Hessian data as a square matrix (unfinished)"
    m_dim = num_atoms(born_file)
    tensor_hess = np.zeros((m_dim,m_dim))
    hess_data = []
    hess = open(hess_file).readlines()
    for i in range(len(hess)):
        hess_split = hess[i].split()
        for j in range(len(hess_split)):
            hess_data.append(float(hess_split[j]))
    start = 0
    end = start + m_dim
    for i in range(m_dim):
        print(start,end)
        print(hess_data[start:end])
        #tensor_hess[i] = hess_data[start:end]
        start = end
        end += m_dim
    return tensor_hess

def parse_hessian(hess_file):
    "Parses Hessian to find index of pairs and diagonal elements"
    hess_data = []
    hess = open(hess_file).readlines()
    for i in range(len(hess)):
        hess_split = hess[i].split()
        for j in range(len(hess_split)):
            if float(hess_split[j])<= 1e-12:
                hess_data.append(0.0)
            else:
                hess_data.append((float(hess_split[j]),i*len(hess_split)+j))
    hess_data = [i for i in hess_data if i != 0.0]
    hess_pairs = []
    hess_diag = []
    flag = True
    while flag == True:
        #print("initial element length ",len(hess_data),"twice pair length + diag length ",2*len(hess_pairs)+len(hess_diag))
        if len(hess_data) == 0:
            flag = False
            break
        for i in range(1,len(hess_data)):
            if abs(hess_data[0][0]-hess_data[i][0]) <= 1e-10:
                mean = (hess_data[0][0]+hess_data[i][0])/2.0
                hess_pairs.append([mean,hess_data[0][1],hess_data[i][1]])
                rem1 = hess_data[0]
                rem2 = hess_data[i]
                hess_data.remove(rem1),hess_data.remove(rem2)
                break
            elif i == len(hess_data)-1:
                hess_diag.append([hess_data[0][0],hess_data[0][1]])
                hess_data.remove(hess_data[0])
    return hess_pairs,hess_diag

def print_indexes(hess_file):
    print("Pair indexes")
    print([i[1:] for i in parse_hessian(hess_file)[0]])
    print("Diag indexes")
    print([i[1] for i in parse_hessian(hess_file)[1]])

print_indexes("Pm.rePhonons.2123667.HESSFREQ.DAT")
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
