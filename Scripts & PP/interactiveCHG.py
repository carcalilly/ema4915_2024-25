#SHELL INTERACTIVE SCRIPT

import numpy  as np
import sys
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import importlib as imp
from   scipy.interpolate import interp1d, make_interp_spline   
from   scipy.integrate   import simpson, trapezoid
from   scipy.signal import savgol_filter
from   scipy.ndimage import gaussian_filter1d

"""
 @author Liping Yu
 heavily edited and repurposed for ema4915 2024/2025
"""

def main():
    global data, a, a1, a2, a3, vol, nat, nx, ny, nz, dx, dy, dz, dvol
    while True:
        chg_file = input("Enter the CHG file name: ").strip() #ask for the input file name
        try:
            with open(chg_file, 'r') as CHG:
                data = CHG.readlines()
            break  # exit loop if successful
        except FileNotFoundError:
            print("File not found. Please try again.")
    print(f"File: {data[0].strip()}") #print the name of the file

    a = float(data[1].split()[0]) #read 3-5 data lines (in data[2-4]), represents lattice vectors in the CHG file 
    a1=a*(np.sum([float(x)**2 for x in data[2].split()]))**0.5 
    a2=a*(np.sum([float(x)**2 for x in data[3].split()]))**0.5
    a3=a*(np.sum([float(x)**2 for x in data[4].split()]))**0.5 #this finds the norm/magnitude of the vector
    vol = a1*a2*a3 #calculates volume (vol) of the unit cell
    print(f"Lattice Parameters: a1={a1:.2f}, a2={a2:.2f}, a3={a3:.2f}")
    print(f"Unit Cell Volume: vol={vol:.2f}")

    nat = np.sum([int(x) for x in data[6].split()]) #line 7 of the CHG file (number of atoms)
    nx,ny,nz = [ int(x) for x in data[nat+9].split() ] #no. of grid points in each direction (x, y, z)
    dx = a1/nx
    dy = a2/ny
    dz = a3/nz #spacing in each direction
    dvol = dx*dy*dz #volume of a grid cell
    print(f"Number of atoms: nat={nat}")
    print(f"Grid Points: nx={nx}, ny={ny}, nz={nz}")
    print(f"Spacing: dx={dx:.2f}, dy={dy:.2f}, dz={dz:.2f}")
    print(f"Grid Cell Volume: dvol={dvol:.2f}")
