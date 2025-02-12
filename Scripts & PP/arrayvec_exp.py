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

for filename in files: #for each file in the folder
    if '-' in filename:
        try:
            at_num = int(filename.split('-')[0]) #pull atomic number from the file name
            z_values.append(at_num) #add atomic number to the list
        except ValueError:
            print(f"Error: {filename} does not follow naming convention.")
        
z_values = np.array(z_values) #convert list to 1D numpy array

print("Atomic Numbers: ", z_values) #print atomic numbers

