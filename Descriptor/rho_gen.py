"""
The arrays and vector representations, experimental script
Use: python array_exp.py directory
"""

#Libraries
import numpy as np
import sys
import os
import subprocess
import matplotlib.pyplot as plt
from   scipy.interpolate import interp1d, make_interp_spline   
from   scipy.integrate   import simpson, trapezoid
from   scipy.signal import savgol_filter
from   scipy.ndimage import gaussian_filter1d

folder_path = sys.argv[1]  # Folder path
output_folder = "outputrrhos"  # output folder
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

#For the input folder that holds the AECCAR files 
folder_path = sys.argv[1] #folder path specified in command line
files = os.listdir(folder_path) #list of files in the folder
elemfolder = [f for f in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, f))]


atomicz = [] #initialize list for atomic numbers
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

atomicm = []
atomic_mass = {
    "H": 1.008, "He": 4.0026, "Li": 7, "Be": 9.012183, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.9984, "Ne": 20.18,
    "Na": 22.98977, "Mg": 24.305, "Al": 26.98154, "Si": 28.085, "P": 30.97376, "S": 32.07, "Cl": 35.45, "Ar": 39.9, "K": 39.0983, "Ca": 40.08,
    "Sc": 44.95591, "Ti": 47.867, "V": 50.9415, "Cr": 51.996, "Mn": 54.93804, "Fe": 55.84, "Co": 58.93319, "Ni": 58.693, "Cu": 63.55, "Zn": 65.4,
    "Ga": 69.723, "Ge": 72.63, "As": 74.92159, "Se": 78.97, "Br": 79.9, "Kr": 83.8, "Rb": 85.468, "Sr": 87.62, "Y": 88.90584, "Zr": 91.22,
    "Nb": 92.90637, "Mo": 95.95, "Tc": 96.90636, "Ru": 101.1, "Rh": 102.9055, "Pd": 106.42, "Ag": 107.868, "Cd": 112.41, "In": 114.818, "Sn": 118.71,
    "Sb": 121.76, "Te": 127.6, "I": 126.9045, "Xe": 131.29, "Cs": 132.9055, "Ba": 137.33, "La": 138.9055, "Ce": 140.116, "Pr": 140.9077, "Nd": 144.24,
    "Pm": 144.9128, "Sm": 150.4, "Eu": 151.964, "Gd": 157.2, "Tb": 158.9254, "Dy": 162.5, "Ho": 164.9303, "Er": 167.26, "Tm": 168.9342, "Yb": 173.05,
    "Lu": 174.9668, "Hf": 178.49, "Ta": 180.9479, "W": 183.84, "Re": 186.207, "Os": 190.2, "Ir": 192.22, "Pt": 195.08, "Au": 196.9666, "Hg": 200.59,
    "Tl": 204.383, "Pb": 207, "Bi": 208.9804, "Po": 208.9824, "At": 209.9872, "Rn": 222.0176, "Fr": 223.0197, "Ra": 226.0254, "Ac": 227.0278,
    "Th": 232.038, "Pa": 231.0359, "U": 238.0289, "Np": 237.0482, "Pu": 244.0642, "Am": 243.0614, "Cm": 247.0704, "Bk": 247.0703, "Cf": 251.0796,
    "Es": 252.083, "Fm": 257.0951}

atomicr = []
atomic_radius_pm = {
    "H": 120, "He": 140, "Li": 182, "Be": 153, "B": 192, "C": 170, "N": 155, "O": 152, "F": 135, "Ne": 154,
    "Na": 227, "Mg": 173, "Al": 184, "Si": 210, "P": 180, "S": 180, "Cl": 175, "Ar": 188, "K": 275, "Ca": 231,
    "Sc": 211, "Ti": 187, "V": 179, "Cr": 189, "Mn": 197, "Fe": 194, "Co": 192, "Ni": 163, "Cu": 140, "Zn": 139,
    "Ga": 187, "Ge": 211, "As": 185, "Se": 190, "Br": 183, "Kr": 202, "Rb": 303, "Sr": 249, "Y": 219, "Zr": 186,
    "Nb": 207, "Mo": 209, "Tc": 209, "Ru": 207, "Rh": 195, "Pd": 202, "Ag": 172, "Cd": 158, "In": 193, "Sn": 217,
    "Sb": 206, "Te": 206, "I": 198, "Xe": 216, "Cs": 343, "Ba": 268, "La": 240, "Ce": 235, "Pr": 239, "Nd": 229,
    "Pm": 236, "Sm": 229, "Eu": 233, "Gd": 237, "Tb": 221, "Dy": 229, "Ho": 216, "Er": 235, "Tm": 227, "Yb": 242,
    "Lu": 221, "Hf": 212, "Ta": 217, "W": 210, "Re": 217, "Os": 216, "Ir": 202, "Pt": 209, "Au": 166, "Hg": 209,
    "Tl": 196, "Pb": 202, "Bi": 207, "Po": 197, "At": 202, "Rn": 220, "Fr": 348, "Ra": 283, "Ac": 260, "Th": 237,
    "Pa": 243, "U": 240, "Np": 221, "Pu": 243, "Am": 244, "Cm": 245, "Bk": 244, "Cf": 245, "Es": 245, "Fm": 257}

ionE = []
ionization_enrg = {
    "H": 13.598, "He": 24.587, "Li": 5.392, "Be": 9.323, "B": 8.298, "C": 11.26, "N": 14.534, "O": 13.618, "F": 17.423, "Ne": 21.565,
    "Na": 5.139, "Mg": 7.646, "Al": 5.986, "Si": 8.152, "P": 10.487, "S": 10.36, "Cl": 12.968, "Ar": 15.76, "K": 4.341, "Ca": 6.113,
    "Sc": 6.561, "Ti": 6.828, "V": 6.746, "Cr": 6.767, "Mn": 7.434, "Fe": 7.902, "Co": 7.881, "Ni": 7.64, "Cu": 7.726, "Zn": 9.394,
    "Ga": 5.999, "Ge": 7.9, "As": 9.815, "Se": 9.752, "Br": 11.814, "Kr": 14, "Rb": 4.177, "Sr": 5.695, "Y": 6.217, "Zr": 6.634,
    "Nb": 6.759, "Mo": 7.092, "Tc": 7.28, "Ru": 7.361, "Rh": 7.459, "Pd": 8.337, "Ag": 7.576, "Cd": 8.994, "In": 5.786, "Sn": 7.344,
    "Sb": 8.64, "Te": 9.01, "I": 10.451, "Xe": 12.13, "Cs": 3.894, "Ba": 5.212, "La": 5.577, "Ce": 5.539, "Pr": 5.464, "Nd": 5.525,
    "Pm": 5.55, "Sm": 5.644, "Eu": 5.67, "Gd": 6.15, "Tb": 5.864, "Dy": 5.939, "Ho": 6.022, "Er": 6.108, "Tm": 6.184, "Yb": 6.254,
    "Lu": 5.426, "Hf": 6.825, "Ta": 7.89, "W": 7.98, "Re": 7.88, "Os": 8.7, "Ir": 9.1, "Pt": 9, "Au": 9.226, "Hg": 10.438,
    "Tl": 6.108, "Pb": 7.417, "Bi": 7.289, "Po": 8.417, "At": 9.5, "Rn": 10.745, "Fr": 3.9, "Ra": 5.279, "Ac": 5.17, "Th": 6.08,
    "Pa": 5.89, "U": 6.194, "Np": 6.266, "Pu": 6.06, "Am": 5.993, "Cm": 6.02, "Bk": 6.23, "Cf": 6.3, "Es": 6.42}


for elemfolder in elemfolder: #for each file in the folder
    CHG = os.path.join(folder_path, elemfolder, "CHG") #path to the file
    try:
        with open(CHG, 'r') as file:
            data = file.readlines() #read lines and store data
        line6 = data[5].strip()
        element = (line6.split('/')[0]) #pull element from line 6
        # print(f"Debug: element variable {element} from {elemfolder}")

        if element in periodic_table:
            at_num = periodic_table[element] #get atomic number from the dictionary
            atomicz.append((at_num)) #add atomic number to the list
            at_mass = atomic_mass[element]
            atomicm.append((at_mass))
            at_rad = atomic_radius_pm[element]  
            atomicr.append((at_rad))
            at_ion = ionization_enrg[element]
            ionE.append((at_ion))

            input_path = CHG # set CHG as input
            output_path = os.path.join(output_folder, f"{element}.txt") # set output to be named by element

            subprocess.run(["python", "calc_rrho_edited.py", input_path, output_path]) # call calc_rrho_edited.py with subprocess
            print(f"Output written to {output_path}")
            
        else: 
            print(f"Error: Cannot read element from {elemfolder}.")
    except ValueError:
            print(f"Error: Problem reading data.")


