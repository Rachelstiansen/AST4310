import numpy as np


def read_NIST_species(input_file):
    """
    Reads data from a text file downloaded from NIST database with
    level energies (4th column and statistical weights (2nd column).
    """
    data = open(input_file, 'r').readlines()[1:]  # skip header
    nlevels = len(data)
    g = np.zeros(nlevels, dtype='i') - 1
    chi = np.zeros(nlevels, dtype='f') - 1
    chi_ion = 0.
    
    for i, line in enumerate(data):
        entries = line.split('\t')
        if entries[1] != '':
            g[i] = int(entries[1])
            chi[i] = float(entries[3].strip('"').strip('[').strip(']'))
        if entries[0].strip('"').lower() == "limit":
            chi_ion = float(entries[3].strip('"').strip('[').strip(']'))
            break
    # clean up missing values
    mask = (g >= 0) & (chi >= 0)
    return g[mask], chi[mask], chi_ion

header = f"#   E (cm^-1)\t\t\tg\tstage" # \t are tabs
output_filename = 'C_I-VI.txt'
with open(output_filename, 'w') as output:
    output.write(header + "\n")
    output.close


for i in range(1, 7): # for every ionization stage 
    filename = f'C_{i}.txt'
    g, chi, chi_ion = read_NIST_species(filename)
    stage = i
'''
    for j in range (len(g)):
        with open(output_filename, 'a') as output:   
            # < is for left alignment, 15 is amount of character space
            append_line = f'\t{chi[j]:<15}\t\t{g[j]}\t{i}\n' # n_stage - 1, due to n_stage = 1 being neutral
            output.write(append_line)
'''       
