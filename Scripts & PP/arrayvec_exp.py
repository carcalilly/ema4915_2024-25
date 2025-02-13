"""
The arrays and vector representations, experimental script
"""

#Libraries
import numpy  as np
import sys
import os



#For the input folder that holds the AECCAR files 
folder_path = sys.argv[1] #folder path specified in command line
files = os.listdir(folder_path) #list of files in the folder


z_values = [] #initialize list for atomic numbers
vasp_file = {} #initialize dictionary for VASP file content 
periodic_table = {"H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
    "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, 
    "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, 
    "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42,
    "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58,
    "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66,
    "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72, "Ta": 73, "W": 74,
    "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82,
    "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90,
    "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98,
    "Es": 99, "Fm": 100} #dictionary for periodic table

#NOTE: something is wrong with the code below, but I have no idea what it is!

for filename in files: #for each file in the folder
    AEC = os.path.join(folder_path, filename) #path to the file
    try:
        open(AEC, 'r')
        element = (AEC[5].split('/')[0]) #pull element from line 6
        
        if element in periodic_table:
            at_num = periodic_table[element] #get atomic number from the dictionary
            z_values.append(at_num) #add atomic number to the list
        else: 
            print(f"Error: Cannot read element from {filename}.")
    except ValueError:
            print(f"Error: Cannot read atomic number from {filename}.")
        
z_values = np.array(z_values) #convert list to 1D numpy array

print("Atomic Numbers: ", z_values) #print atomic numbers


"""
NOTE: code for extracting atomic numbers from file name
for filename in files: #for each file in the folder
    if '-' in filename:
        try:
            at_num = int(filename.split('-')[0]) #pull atomic number from the file name
            z_values.append(at_num) #add atomic number to the list
        except ValueError:
            print(f"Error: {filename} does not follow naming convention.")
        
z_values = np.array(z_values) #convert list to 1D numpy array

print("Atomic Numbers: ", z_values) #print atomic numbers
"""


