import numpy as np
import sys
import csv
import os

# python Descriptorinator.py OUTPUTFILENAME

input_folder = "outputrrhos"  # directory
output_file = sys.argv[1]  # output file
while not os.path.exists(input_folder):
    print(f"Error: Folder '{input_folder}' not found.")
    input_folder = input("Please enter a valid folder path containing calc_rrho outputs: ")

fout = open(output_file, 'w')
Desc = []

# List all files in the output folder and sort them
entries = sorted(os.listdir(input_folder))

for entry in entries:
    file_path = os.path.join(input_folder, entry)
    ext = os.path.splitext(entry)[-1].lower()
    
    if ext == '.txt':  # process only .txt files
        try:
            with open(file_path, 'r') as file:
                rawdata = file.read().split()[12::4]  # every 4th string starting from the 13th
            
            idata = [float(rawdata[i]) for i in range(300)]  # First 300 points
            jdata = []
            
            for j in range(50):  # 50 points in the descriptor array
                bob = 0.01 * sum(idata[j * 6:(j + 1) * 6])  # sum 6 points, multiply by dr=0.01
                fout.write(f"{bob:15.8E}\t")
                jdata.append(bob)
            
            Desc.append(jdata)
            fout.write("\n")

        except (IndexError, ValueError, FileNotFoundError) as e:
            print(f"Error processing file {entry}: {e}")

fout.close()

# convert descriptor list into a 2D NumPy array
if Desc:
    TwoDscriptor = np.reshape(np.array(Desc), (len(Desc), 50))
    print("Descriptor array successfully created.")
else:
    print("No valid data processed.")
