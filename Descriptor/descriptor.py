import numpy as np
import sys
import csv
import os
#python3 Descriptorinator.py OUTPUT
#Make sure all elements have same file extension and nothing else in directory share that extension
#Each element has its own row, Ascending order; 1->92 from top-> bottom
#Read element's file -> create its array -> add to the descriptor

fout = open(sys.argv[1], 'w')
Desc=[]
with os.scandir() as entries:               #ALL files in current directory
    sentries = sorted(entries, key=lambda e:e.name) #Sort alphanumerically (os.scandir is random)
    for entry in sentries:                  #Repeat following for every file in directory
        ext = os.path.splitext(entry)[-1].lower()
        if ext == '.txt':                   #ONLY .txt files allowed
            rawdata = open(entry, 'r').read().split()[12::4]
            #rawdata starts at 13th string of the input file, and every 4th string beyond; i.e. this is only for calc_rrho outputs
            idata = [] #Empty arrays for later
            jdata = [] #>:(
            
            for i in range(300): #First 300 points; r increasing by 0.01 per point (First 3 ANGSTROM)
                idata.append(float(rawdata[i])) #Convert ith STRING of rawdata to ith FLOAT of idata

            for j in range(50): #50 Points in the descriptor array; add 300/50 = 6 "real" points per descriptor point
                bob=0.01*float(idata[j*6] + idata[j*6+1] + idata[j*6+2] + idata[j*6+3] + idata[j*6+4] + idata[j*6+5])
                    #Multiply e/A by dr = 0.01A to get e             ^There must be a better way to do this
                fout.write("%15.8E\t" %(bob)) #Data to file
                jdata.append(bob) #50 bobs form the list jdata
            Desc.append(jdata) #Add 1 element to descriptor (still a list!)
            fout.write("\n") #New line for next element
fout.close()

TwoDscriptor = np.reshape(np.array(Desc),(92,50)) #Turn list into proper 2D array
#Maybe we just continue from here instead of making another script that reads this one's output
