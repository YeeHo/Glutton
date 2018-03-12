# Glutton
Part 1 of 3. Use Glutton script to generate structure ensemble based on chemical shifts

Step1. Put database files (LEVEL*.dat), Glutton_R1.0.py and input.txt files in the same folder.

Step2. create subfolders named "cs" and "PDBOUT" within this folder.

Step3. Put input chemical shift file (NMRSTAR 3.0 format) in the subfolder "cs"

Step4. change the input.txt to set parameters to determine the candidate pool

        bmr15398.str   # name of the input chemical shift file 
        2              # LEVEL of the database, 1 - high resolution; 2 - medium resolution; 3 - low resolution
        10             # number of structures to be generated
        0.5            # the selected width of the chemical shift distribution used to derive statistical distributions 
        False          # use log of the densisty distributions of Phi and Psi 
                       # angles as probability function to generate angles

        NOTE: if there are insufficient data to generate distributions, the choice can be 
        [1]. use low resolution database - LEVEL3;            
        [2]. increase the selected width of the chemical shift distribution. 
 
Step5. run script as following:

    python Glutton_R1.0.py
 
Part 2 of 3. Libraries needed to run Glutton script

The Glutton script has been tested using Anaconda python distribution version 3.6.

The Anaconda python package can be downloaded for free at: https://www.anaconda.com/download

Additional libraries needed to run Glutton script:

[1]. Install Biopython

        conda install biopython (or use "pip install biopython")
  
  Biopython webpage: http://biopython.org/
  
[2]. Install PeptideBuilder

        pip install PeptideBuilder

  PeptideBuilder webpage: https://github.com/mtien/PeptideBuilder
  
  Reference:
  M. Z. Tien, D. K. Sydykova, A. G. Meyer, C. O. Wilke (2013). PeptideBuilder:
  A simple Python library to generate model peptides. PeerJ 1:e80.

Part 3 of 3. Protein Chemical Shift - Structure Database 

All database files 

LEVEL1_ALL_20180116.dat

LEVEL2_ALL_20180116.dat

LEVEL3_ALL_20180116.dat

are in csv format. 

Each row represents a single residue, and the meaning of each column is as follows:

column 1: residue type (one letter code)

column 2: the chemical shift value of the H on the backbone N

column 3: the chemical shift value of the backbone N

column 4: the chemical shift value of the H on the backbone C_Alpha

column 5: the chemical shift value of the backbone C_Alpha

column 6: the chemical shift value of the backbone C'

column 7: Secondary structure type of this residue

column 7: Phi angle of this residue

column 7: Psi angle of this residue

NOTE: if the chemical shift value of a atom is not avaliable, Glutton sets it as 9999.0. 
