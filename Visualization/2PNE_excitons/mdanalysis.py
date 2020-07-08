import MDAnalysis as md
import numpy as np

#Load the pdb file for 2PNE

u = md.Universe("2pne.pdb")

#Create an emtpy position array for the carbons
a=np.zeros((84,3))
Carbon = u.select_atoms("name C")
#Specify the sidechains
side = u.select_atoms("resname ASN and name CG")
compleet=[Carbon.positions,side.positions]
a[0:81,:]=Carbon.positions
a[80:84,:]=side.positions
np.savetxt("position.txt",a)
