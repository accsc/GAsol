#!/bin/bash

# Run GAsol with the 3D-RISM density in complex.O.1.dx 
# and in a restriced area given by a ligand (lig.pdb) and radius of 10.0 A
# and a ratio of densitiy/radius for each water molecule of more than or equal to 0.25
# to filter out really weak waters
time ../../gasol --dx complex.O.1.dx --lig lig.pdb -r 8.0 > solution_gasol.pdb
# The program is able to identify all conserved water molecules in BRD4 (5I80_A_waters_ref.pdb)
