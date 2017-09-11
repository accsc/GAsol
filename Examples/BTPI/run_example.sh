#!/bin/bash

# Run GAsol with the 3D-RISM density in complex.O.1.dx 
# and in a restriced area given by a center (31.2,6.8,1.7) and radius of 8.0 A
time ../gasol --dx complex.O.1.dx -x 31.242  -y 6.845  -z 1.776 -r 8.0 > solution_gasol.pdb
# The program is able to identify all three water molecules in the cavitiy
