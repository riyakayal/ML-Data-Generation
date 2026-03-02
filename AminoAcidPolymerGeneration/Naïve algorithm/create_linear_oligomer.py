# Author: Riya Kayal
# Created: 10/10/2025

# Generate an oligomer by adding n number of monomers on a stright line.
# NOTE: This is not a blackbox program. Find your reference atom (ATOM1-1) 
# in first molecule (mol1) and the distance 
# of it from the same reference atom (ATOM1-2) in the second molecule (mol2).
# We take glycine as an example, to add n number of glycines to it.
# We place it as *HOOC-CH2-NH2 HOOC-CH2-NH2*, keeping the NH2 group of the first molecule 
# in a suitable distance to make H bond formation possible between the carboxylic O atom
# of the second molecule with the H from the first molecule's NH2 group.
# Find the O (ATOM1-1) in the C=O of COOH, here it is the 4th one (m = mol1[3]). 
# This O atom in secnd glycine should be close to N (ATOM2-1) of first glycine. 
# get this O's coordinate by drawing it about 2.5 A away from N of first glycine in Avogadro.

mol2_O = "O       -0.7049115307     -1.1265436717      1.8039367523"
mol2_O = mol2_O.split()

# Now we know O's (from C=O) coordinates in both glycines.
# Translate first glycine by the distance between these two Os.
# Repeat the process to get the oligomer of desired size.
'''
mol1 = ["N    1.3476865998     -1.6766316554      1.1486751945",             
"C        2.4805893087     -0.7619959339      1.1783488893",                 
"C        3.6383064077     -1.4243062858      0.4375805008",                 
"O        4.7865995315     -1.0749215957      0.5402265471",                 
"O        3.2608914825     -2.4292203189     -0.3357786539",                 
"H        0.4894402933     -1.2064395997      1.4097403346",                 
"H        1.5092168303     -2.4247806682      1.8221827028",                 
"H        2.8241714092     -0.4747813693      2.1802151016",                 
"H        2.2274992088      0.1527784179      0.6325934985",                 
"H        2.2899993736     -2.5464142536     -0.1975819244"]
'''

# read the coordinate file
# get mol1 in the above format
def get_translation_vector(xyzfile):
    molecule1 = []
    f = open(xyzfile, 'r') 
    for line in f.readlines()[2:]:
        line = line.strip() 
        molecule1.append(line)
    return molecule1

# create oligomer
def get_oligomers(n_oligomers, monomer_name, xyzfile):
    
    print("="*50, "\n Preparing oligomer coordinates for %s..." % monomer_name)
    print("="*50, "\n")
    mol1 = get_translation_vector(xyzfile)
    # We choose the first O atom as the reference atom (ATOM1)
    m = mol1[3].split()
    mx, my, mz = float(mol2_O[1]), float(mol2_O[2]), float(mol2_O[3])
    # translation vectors along x, y & z axes
    d_x, d_y, d_z = float(m[1]) - mx, float(m[2]) - my, float(m[3]) - mz; #print(d_x)
    dict1 = {2: "2nd", 3: "3rd"}
    
    for times in range(1, n_oligomers):
        print("-"*50,"\n Coordinates of", dict1.get(1 + times, "%2dth"%(1 + times)), monomer_name)
        print("-"*50)
        for atom in mol1:
            element, x, y, z = atom.split()
            [x,y,z] = list(map(float, [x,y,z]))
            #translate 1st monomer's coordinates by (i -1) * d_j where i indicates i-th oligomer, d_j is displacement with j = {x,y,z}
            x -= times * d_x
            y -= times * d_y
            z -= times * d_z
            print(element, "  %10.10f  %10.10f  %10.10f" % (x, y, z))
        print("\n")
        
# run         
if __name__ == "__main__":
    print("Enter name of the .xyz file without the .xyz suffix:")
    name = input()
    print("Enter oligomer length as natural number:")
    chain_length = int(input())
    print("\n")
    get_oligomers(chain_length, name, name + ".xyz")
