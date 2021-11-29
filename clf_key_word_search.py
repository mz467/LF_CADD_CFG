# Coupled LAMMPS and FEniCS
# To read in the standard input 
# Written by Wenjia G. on 08/31/2018

import sys

def key_search(string, input):
    string = string.lower()
    for i in range(len(input)):
        line = input[i].strip().lower()
        if string in line:
            line = line.split()
            var = line[0]
            break
    if 'var' not in locals():
        print("ERROR: %s cannot be found in the input" % (string))
        exit()
    return var



