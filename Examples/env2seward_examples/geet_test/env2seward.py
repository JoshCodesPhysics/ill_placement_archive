#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 15:33:04 2020

@author: joshhorswill10
"""
import sys
import json

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
    for i in range(len(pseudo_data)):
        if pseudo_data[i].split() == ['En','symetrie'] or pseudo_data[i].split() == ['en','symetrie']:
            pseudo_data = pseudo_data[i+1:]
            break
    for i in range(len(frag_data)):
        element = frag_data[i].split()[0]
        for j in range(len(element)):
            if element[j].isdigit():
                atom = element[0:j]
                if atom not in frag_list:
                    frag_list.append(atom)
                break
            elif element[j] == '_':
                atom = element[0:j]
                if atom not in frag_list:
                    frag_list.append(atom)
    for i in range(len(pseudo_data)):
        element = pseudo_data[i].split()[0]
        for j in range(len(element)):
            if element[j].isdigit():
                atom = element[0:j+1]
                if atom not in pseudo_list:
                    pseudo_list.append(atom)
                break
            elif element[j] == '_':
                atom = element[0:j]
                if atom not in pseudo_list:
                    pseudo_list.append(atom)
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
       sym_bool = False
       frag = open(sew0_file).readlines()
       #Create list of frag data
       for i in range(len(frag)):
           if frag[i].split()[0] == 'SYMMETRY':
               sym_type = frag[i+1].split()[0]
               sym_bool = True
           elif frag[i].split()[0][:8] == "Fragment":
               start = i+1
           elif frag[i].split()[0][:7] == "Pseudos":
               end = i
               break
       frag_sew0 = frag[start:end]
       for i in range(len(frag_sew0)):
           line_list = frag_sew0[i].split()
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
       if sym_bool == True:
           sym_list = [' SYMMETRY\n',' %s\n'%sym_type,'\n']
           sew_in.writelines(sym_list)
       sew_in.write('*** Fragment ***************************************************\n')
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
                for i in range(len(psd_whole)):
                    if psd_whole[i].split() != []:    
                        psd_split = psd_whole[i].split()
                        if psd_split == ['En','symetrie'] or psd_split == ['en','symetrie']:
                            new_start = i+1
                            psd_sym = psd_whole[new_start:]
                            psd_whole = psd_sym
                            break
                psd = []
                for i in range(len(psd_whole)):
                    if psd_whole[i].split()[0] == 'Atom':
                        pass
                    elif psd_whole[i].split()[0][1].isdigit() == False:
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
                    #print("element",elem,"coords", coords,"psd_coords", psd_coords)
                    separation = ((psd_coords[0]-coords[0])**2+(psd_coords[1]-coords[1])**2+(psd_coords[2]-coords[2])**2)**0.5
                    #print("Separation ", separation)
                    if separation <= 1e-1 and elem == psd[i].split()[0][0:2]:
                        if '_' in psd[i].split()[0]:
                            atom = psd[i].split()[0][0:psd[i].split()[0].index('_')]
                        else:
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

def atom_print(atoms):
    print("Fragment atoms: ",atoms[0], "TIP atoms: ", atoms[1])

def finalprompt():
    write_prompts = ask_user('Would you like to save these prompts to an input file?')
    title = str(input('Please enter TITLE line (Maximum 80 characters):   '))
    prefix = str(input('Please enter files PREFIX (must match input file prefixes):   '))
    sew0name = "%s.env.sew0"%prefix
    psdname = "%s.env.psd"%prefix
    filename = "%s.sew.in"%prefix
    atoms = atom_find(sew0name,psdname)
    atom_print(atoms)
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
    if write_prompts == True:
        input_file = open('%s.in'%prefix,'w')
        input_list = ['filename = %s\n'%filename,'title = %s\n'%title,\
                'sew0_file = %s\n'%sew0name,'psd_file = %s\n'%psdname, 'lib_frag = %s\n'%str(lib_frag),\
                'lib_pseudo = %s'%str(lib_pseudo)]
        input_file.writelines(input_list)
        input_file.close()
        print("Input file generated")
    print("Inputs complete")
    finalwrite(filename,title,sew0name,psdname,lib_frag,lib_pseudo)

def fileinput(input_file):
    file = open(input_file).readlines()
    for i in range(len(file)):
        file_line = file[i].split()
        if file_line[0] == "filename":
            filename = file_line[2]
        elif file_line[0] == "title":
            title = file_line[2]
        elif file_line[0] == "sew0_file":
            sew0_file = file_line[2]
        elif file_line[0] == "psd_file":
            psd_file = file_line[2]
        elif file_line[0] == "lib_frag":
            lib_frag = "".join(file_line[2:])
            lib_frag = lib_frag.replace("'",'"')
            lib_frag = json.loads(lib_frag)
        elif file_line[0] == "lib_pseudo":
            lib_pseudo = "".join(file_line[2:])
            lib_pseudo = lib_pseudo.replace("'",'"')
            lib_pseudo = json.loads(lib_pseudo)
    atoms = atom_find(sew0_file,psd_file)
    atom_print(atoms)
    finalwrite(filename,title,sew0_file,psd_file,lib_frag,lib_pseudo)

def input_jupyter_or_prompt():
   if len(sys.argv) >= 2:
       fileinput(sys.argv[1])
   else:
       if ask_user("Are you running this on the Jupyter notebook?") == False:
          finalprompt()
input_jupyter_or_prompt()
