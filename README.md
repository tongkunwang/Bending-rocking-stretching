# Bending-rocking-stretching
This stupid code parses vibrational modes of O atoms in zeolite frameworks into bending, rocking and stretching by a very inefficient way.
It is (poorly) designed for the format of VASP (a periodic DFT simulation package) and VESTA (an visualization software) to study vibrations on O atoms in zeolite frameworks.

This version can only handle pure-silica zeolite frameworks (only O and Si atoms), without any other species like OSDA.

To use this code, we need 3 files in the same folder:
1) This code
2) The POSCAR file of the system of interest
3) A file containing the "VECTR", which contains vibrational vectors in vesta format. (In this code the file is invoked also as VECTR, but you can change that)

Examples of a POSCAR file, a file named VECTR and the original vesta file that contains the VECTR part are in the list of files here.

The original vesta file was generated by extracting vibrational modes from an OUTCAR file. There are some existing codes that can do this, like https://github.com/Stanford-MCTG/VASP-plot-modes

Headers of the code need to be changed specifically for structural details of the zeolite system you want to calculate
1) Lattice parameters: axis_a, axis_b, axis_c
2) Axis angles: alpha, beta, gamma
3) Numbers of atoms: o_number, si_number, atom_number

All parameters and files names are fixed in the source code so there is nothing we can put in through command line...
Please change the code to change parameters and file names and let's blame it on my trivial programming skill and bad habbit of alway writting new codes on old codes.(O_O)

Finally we will get data in columns of index of atoms, Si-O-Si bond angles, components for decomposed bending and rocking modes, and modules of vibrational vectors.
