#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 15:33:04 2020

@author: joshhorswill10
"""
import sys

def opening_lines(new_file,title_input):
    """This function writes initial opening commands before the basis section and
    returns the new file with title from the function input"""
    #Editing new titles to name new file and write third line
    sew_in = open(new_file,'w')
    opening_list = ["&SEWARD\n","Title\n"," "+title_input+"\n","\n","Expert\n","Verbose\n","\n"]
    sew_in.writelines(opening_list)
    sew_in.close()

def atom_find(sew0_file,psd_file):
    frag_list = []
    pseudo_list = []
    for i in range(len(open(sew0_file).readlines())):
        if open(sew0_file).readlines()[i].split()[0] == "Fragment":
            start = i+1
        elif open(sew0_file).readlines()[i].split()[0] == "Pseudos":
            end = i
            break
    frag_data = open(sew0_file).readlines()[start:end]
    pseudo_data = open(psd_file).readlines()[1:]
    for i in range(len(frag_data)):
        element = frag_data[i].split()[0]
        for j in range(len(element)):
            if element[j].isdigit():
                atom = element[0:j]
                if atom not in frag_list:
                    frag_list.append(atom)
                break
    for i in range(len(pseudo_data)):
        element = pseudo_data[i].split()[0]
        for j in range(len(element)):
            if element[j].isdigit():
                atom = element[0:j+1]
                if atom not in pseudo_list:
                    pseudo_list.append(atom)
                break
            elif j == len(element)-1:
                atom = element
                if atom not in pseudo_list:
                    pseudo_list.append(atom)
    frag_list.sort(),pseudo_list.sort()
    return frag_list,pseudo_list

def frag_basis(new_file,sew0_file,psd_file,lib):
    "This formats the fragment basis into groups, readable by the next piece of software"
    with open(new_file,'a') as sew_in:
       frag_dict = {}
       atom_list = atom_find(sew0_file,psd_file)[0]
       atom_check = []
       new_atom = False
       frag = open(sew0_file).readlines()
       #Create list of pseudo data
       for i in range(len(frag)):
           if frag[i].split()[0] == "Pseudos":
               end = i
               break
       frag_sew0 = frag[:end]
       for i in range(len(frag_sew0)):
           line_list = frag_sew0[i].split()
           if line_list[0] == "Fragment":
               pass
           else:
               atom_item = line_list[0]
               #print("atom_item: ",atom_item)
               for j in range(len(atom_item)):
                   if atom_item[j].isdigit():
                       atom = atom_item[:j]
                       #print("atom: ",atom)
                       if atom not in atom_check:
                           new_atom = True
                           #print("atom boolean: ",new_atom)
                           atom_check.append(atom)
                       else:
                           new_atom = False
                           #print("atom boolean: ",new_atom)
                       break
               if new_atom == True:
                   frag_dict[atom] = []
                   data_list = frag_sew0[i]
                   frag_dict[atom].append(data_list)
               else:
                   data_list = frag_sew0[i]
                   frag_dict[atom].append(data_list)
       for i in range(len(atom_list)):
           sew_in.write("Basis set\n")
           if lib[atom_list[i]]["loc"] == "":
               sew_in.write(" %s\n"%lib[atom_list[i]]["key"])
           else:
               sew_in.write(" %s    / %s\n"%(lib[atom_list[i]]["key"],lib[atom_list[i]]["loc"]))
           sew_in.write("  spherical\n")
           sew_in.writelines(frag_dict[atom_list[i]])
           sew_in.write("End of basis\n")
           sew_in.write("****\n")

def pseudos(new_file,sew0_file,psd_file,lib):
    "Format of ion pseudopotential"
    "Enter key string as usual except put $ in as follows:"
    "Gd.ECP.Marie.0s.0s.0e-$-GdMn2O5. or similarly Mn.ECP.Marie.0s.0s.0e-$-GdMn2O5."
    "Code will replace this with the corresponding atom number, do not enter / PSEUDO part"
    with open(new_file,'a') as sew_in:
        #Dictionary to store positional data
        pseudo_dict = {}
        atom_list = atom_find(sew0_file,psd_file)[1]
        atom_check = []
        sew0 = open(sew0_file).readlines()
        new_atom = False
        #Create list of pseudo data
        for i in range(len(sew0)):
            if sew0[i].split()[0] == "Pseudos":
                start = i
            elif sew0[i].split()[0] == "XFIEld":
                end = i
                break
        pseudo_sew0 = sew0[start:end]
        #Scan the list
        for j in range(len(pseudo_sew0)):
            #Ignore first line
            if pseudo_sew0[j].split()[0] == "Pseudos":
                pass
            #If the atom is one character, it's a lot more simple
            elif pseudo_sew0[j].split()[0][1].isdigit():
                atom = pseudo_sew0[j].split()[0][0:2]
                #Is it new?
                if atom not in atom_check:
                    new_atom = True
                    atom_check.append(atom)
                else:
                    new_atom = False
                #Append info to dictionary
                if new_atom == True:
                    pseudo_dict[atom] = []
                    pseudo_dict[atom].append(pseudo_sew0[j])
                else:
                    pseudo_dict[atom].append(pseudo_sew0[j])
            #It has two characters, oh dear
            else:
                psd_whole = open(psd_file).readlines()
                psd = []
                for i in range(1,len(psd_whole)):
                    if psd_whole[i].split()[0][1].isdigit() == False:
                        psd.append(psd_whole[i])
                elem = pseudo_sew0[j].split()[0][0:2]
                coords = []
                for i in range(1,4):
                    coords.append(float(pseudo_sew0[j].split()[i]))
                #scan psd for corresponding coords
                for i in range(len(psd)):
                    psd_coords = []
                    for k in range(1,4):
                        psd_coords.append(float(psd[i].split()[k]))
                    separation = ((psd_coords[0]-coords[0])**2+(psd_coords[1]-coords[1])**2+(psd_coords[2]-coords[2])**2)**0.5
                    if separation <= 1e-1 and elem == psd[i].split()[0][0:2]:
                        atom = psd[i].split()[0][0:3]
                        break
                #Is it new?
                if atom not in atom_check:
                    new_atom = True
                    atom_check.append(atom)
                else:
                    new_atom = False
                #Append info to dictionary
                if new_atom == True:
                    pseudo_dict[atom] = []
                    pseudo_dict[atom].append(pseudo_sew0[j])
                else:
                    pseudo_dict[atom].append(pseudo_sew0[j])
        for i in range(len(atom_list)):
            for k in range(len(pseudo_dict[atom_list[i]])):
                if pseudo_dict[atom_list[i]][k].split()[0][1].isdigit():
                    if pseudo_dict[atom_list[i]][k].split()[0][0:2] != atom_list[i]:
                        sys.exit("Atoms have been grouped incorrectly. Check that .sew0 and .psd match properly.")
                else:
                    if pseudo_dict[atom_list[i]][k].split()[0][:2] != atom_list[i][:2]:
                        sys.exit("Atoms have been grouped incorrectly. Check that .sew0 and .psd match properly.")
        sew_in.write("*** Pseudos ***************************************************\n")
        for i in range(len(atom_list)):
            sew_in.write("Basis Set\n")
            if lib[atom_list[i]]["loc"] == "":
                sew_in.write(" %s\n"%lib[atom_list[i]]["key"])
            else:
                sew_in.write(" %s    / %s\n"%(lib[atom_list[i]]["key"],lib[atom_list[i]]["loc"]))
            sew_in.write("  pseudocharge\n")
            data_list = pseudo_dict[atom_list[i]]
            sew_in.writelines(data_list)
            sew_in.write("End of Basis\n")
            sew_in.write("*\n")
        sew_in.write("***************************************************************\n")

def xfield(new_file,sew0_file):
    "A direct copy and paste of the Madelung potential from the .env.sew0 file after scanning"
    with open(new_file,'a') as sew_in:
        num_lines = len(open(sew0_file).readlines())
        for i in range(num_lines):
            if open(sew0_file).readlines()[i] == "XFIEld\n":
                start = i
                break
            elif i == num_lines:
                print("Xfield not found")
        xfield = open(sew0_file).readlines()[start:num_lines]
        xfield.append("End of input\n")
        xfield.append("\n")
        sew_in.writelines(xfield)

def finalwrite(new_file,title_input,sew0_file,psd_file,lib_frag,lib_pseud):
    opening_lines(new_file,title_input)
    frag_basis(new_file,sew0_file,psd_file,lib_frag)
    pseudos(new_file, sew0_file, psd_file,lib_pseud)
    xfield(new_file,sew0_file)
    print("File has been created")

def ask_user(question):
    check = str(input("%s (Y/N): "%question)).lower().strip()
    try:
        if check[0] == 'y':
            return True
        elif check[0] == 'n':
            return False
        else:
            print('Invalid Input')
            return ask_user(question)
    except Exception as error:
        print("Please enter valid inputs")
        print(error)
        return ask_user(question)

def finalprompt():
    title = str(input('Please enter TITLE line (Maximum 80 characters):   '))
    prefix = str(input('Please enter files PREFIX (must match input file prefixes):   '))
    sew0name = "%s.env.sew0"%prefix
    psdname = "%s.env.psd"%prefix
    filename = "%s.sew.in"%prefix
    atoms = atom_find(sew0name,psdname)
    lib_frag = {key: {} for key in atoms[0]}
    lib_pseudo = {key: {} for key in atoms[1]}
    print("")
    print("FRAGMENT ATOMS BASIS:")
    for i in range(len(atoms[0])):
        input_atom = atoms[0][i]
        lib_frag[input_atom]["key"] = input("Please enter the basis set for %s:   "%input_atom)
        lib_frag[input_atom]["loc"] = input("LIBRARY LOCATION (leave blank if default):  ")
    print("")
    print("TIPS ATOMS BASIS:")
    for i in range(len(atoms[1])):
        input_atom = atoms[1][i]
        lib_pseudo[input_atom]["key"] = input("Please enter the TIPS for %s:   "%input_atom)
        lib_pseudo[input_atom]["loc"] = input("LIBRARY LOCATION (leave blank if default):  ")
    print("Inputs complete")
    finalwrite(filename,title,sew0name,psdname,lib_frag,lib_pseudo)

"For testing:"
"""
opening_lines("yes","GdMn2O5_J1.sew.in 2")
lib1 = {'Mn':{'bool': False,'loc': '','key':'Mn.ano-rcc.Roos.21s15p10d6f4g2h.6s4p3d1f0g.'},'O':{'bool': False,'loc': '','key':'O.ano-rcc.Roos.14s9p4d3f2g.4s3p1d0f'}}
frag_basis("yes","GdMn2O5_J1.env.sew0","GdMn2O5_J1.env.psd",lib1)
lib2 = {'Gd1':{'bool': True,'loc': 'PSEUDO','key':'Gd.ECP.Marie.0s.0s.0e-Gd1-GdMn2O5.'},'Gd2':{'bool': True,'loc': 'PSEUDO','key':'Gd.ECP.Marie.0s.0s.0e-Gd2-GdMn2O5.'},'Mn1':{'bool':True,'loc':'PSEUDO','key':'Mn.ECP.Marie.0s.0s.0e-Mn1-GdMn2O5.'},'Mn2':{'bool':True,'loc':'PSEUDO','key':'Mn.ECP.Marie.0s.0s.0e-Mn2-GdMn2O5.'},'O1':{'bool':True,'loc':'PSEUDO','key':'O.ECP.Marie.0s.0s.0e-O1-GdMn2O5.'},'O2':{'bool':True,'loc':'PSEUDO','key':'O.ECP.Marie.0s.0s.0e-O2-GdMn2O5.'},'O3':{'bool':True,'loc':'PSEUDO','key':'O.ECP.Marie.0s.0s.0e-O3-GdMn2O5.'},'O4':{'bool':True,'loc':'PSEUDO','key':'O.ECP.Marie.0s.0s.0e-O4-GdMn2O5.'}}
pseudos("yes","GdMn2O5_J1.env.sew0","GdMn2O5_J1.env.psd",lib2)
xfield("yes","GdMn2O5_J1.env.sew0")
"""

def jupyter_or_prompt():
   if ask_user("Are you running this on the Jupyter notebook?") == False:
      finalprompt()
   else:
      pass
jupyter_or_prompt()
