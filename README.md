All the input files and scripts used to create the data and analyze them for my thesis titled:

"Examining the Internal Structure of Protein Aggregates using Two-Dimeansional Infrared Spectroscopy"

All MD analysis scripts have been adapted from the scripts used with the "Lyzosime in water" GROMACS tutorial [1]

The InpTra files, located in the AmideMapping directory, are used to translate the energy files created by the Amide-I mapping software.[2] This means that time intervals as well as the output format can be changed as specified. This is used for the isotopic labeling and an example input file, showing the necessary formating, has been given. After translation the binary energy and dipole files can be used within NISE.

The exciton analysis script, used on the 2PNE protein, was based on previous work by T.L.C. Jansen [3] with the help of the MD analysis Software package [4]

The 6Y1A isolabel graphing script runs over a set of directories that contain the data from NISE for the isotopically labeled samples. An average Hamiltonian created with the anaylisis function of NISE is necessary, as wel as the "TD_Absorption" and "Absorption.dat" files which correspond to a couplingcut set to zero (Fully labeled), and "TD_Absorption_2" and "Absorption_2.dat" files which are the result of the couplingcut set to 100 (Diluted sample).



[1] Justin A. Lemkul, "Gromacs Tutorial" ,http://www.mdtutorials.com/gmx/lysozyme/index.html

[2] T.L.C. Jansen, "NISE_2017" https://github.com/GHlacour/NISE_2017

[3] https://github.com/lacourjansenlab/Hamiltonians/tree/master/Spectra

[4] https://www.mdanalysis.org/
