"""
The arrays and vector representations, experimental script
Use: python array_exp.py directory
"""

#Libraries
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from   scipy.interpolate import interp1d, make_interp_spline   
from   scipy.integrate   import simpson, trapezoid
from   scipy.signal import savgol_filter
from   scipy.ndimage import gaussian_filter1d



#For the input folder that holds the AECCAR files 
folder_path = sys.argv[1] #folder path specified in command line
files = os.listdir(folder_path) #list of files in the folder
elemfolder = [f for f in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, f))]


zvschg = [] #initialize list for atomic numbers
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

def calcrrho(data):
    
    a = float(data[1].split()[0])
    a1=a*(np.sum([float(x)**2 for x in data[2].split()]))**0.5 
    a2=a*(np.sum([float(x)**2 for x in data[3].split()]))**0.5
    a3=a*(np.sum([float(x)**2 for x in data[4].split()]))**0.5
    vol = a1*a2*a3

    nat      = np.sum([int(x) for x in data[6].split()])
    nx,ny,nz = [ int(x) for x in data[nat+9].split() ]
    dx = a1/nx
    dy = a2/ny
    dz = a3/nz
    dvol = dx*dy*dz

    ncolums    = len(data[nat+11].split())
    grid_lines = int(nx*ny*nz/ncolums)
    if float(nx*ny*nz)/ncolums > grid_lines :
        grid_lines += 1

    nat_lines = int(nat/ncolums)
    if float(nat)/ncolums > nat_lines : 
        nat_lines  += 1

    data_1D = []
    for l in range(nat+10, nat+10+grid_lines) : 
        for x in data[l].split() :              
            data_1D.append(float(x))
    data_3D = np.reshape(np.array(data_1D),(nx,ny,nz),order='F')/vol

    result = np.sum(data_3D)*dvol
    return result

for elemfolder in elemfolder: #for each file in the folder
    CHG = os.path.join(folder_path, elemfolder, "CHG") #path to the file
    try:
        with open(CHG, 'r') as file:
            data = file.readlines() #read lines and store data
        line6 = data[5].strip()
        element = (line6.split('/')[0]) #pull element from line 6
        # print(f"Debug: element variable {element} from {elemfolder}")

        tot_charge = calcrrho(data)

        if element in periodic_table:
            at_num = periodic_table[element] #get atomic number from the dictionary
            zvschg.append((at_num, tot_charge)) #add atomic number to the list
        else: 
            print(f"Error: Cannot read element from {elemfolder}.")
    except ValueError:
            print(f"Error: Cannot read atomic number from {elemfolder}.")
        
z_totchg_arr = np.array(zvschg) #convert list to 2D numpy array
z_totchg_arr = z_totchg_arr[z_totchg_arr[:,0].argsort()]
print(f"Atomic Z vs. Charge \n", z_totchg_arr)


