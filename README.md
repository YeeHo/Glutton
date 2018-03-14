# Glutton
Please take 5 minutes to read this brief readme file before running the Glutton script. 

Part 1 of 4. Introduction of Glutton

        The Glutton database that we have developed is specifically designed to fill a large void
        currently existing in the bioinformatics tools for the structural analysis of intrinsically
        disordered proteins (IDPs) based on experimental information, namely NMR. Our Glutton database 
        contains a total of 5,270 proteins. Glutton database has a hierarchical system of three tiers
        (levels of resolution) that permit customization of the database for specific applications: 
        for structural ensembles of IDPs, or for protein structure determination and prediction of 
        well folded structures. And, most importantly, Glutton changes the structure classification 
        paradigm from the conventional one (one set of chemical shifts to one structure) to one that 
        interprets chemical shifts as ensemble averages, and thus calculates conformational distributions 
        from chemical shift datasets. These distributions permit the quick determination of maximally broad 
        structural ensembles consistent with the data, and thus optimal sampling of the conformational space 
        available to the IDP. These ensembles provide a realistic description of the conformational 
        heterogeneity and structural propensities of IDPs, but they can also be used as starting point 
        for the refinement of protein structure determination (adding additional experimental information)
        or protein prediction (incorporating a molecular force-field).

Part 2 of 4. Use Glutton script to generate structure ensembles based on chemical shifts

        Step1. Put database files (LEVEL*.dat), Glutton_R1.0.py and input.txt files in the same folder.

        Step2. If there are no subfolders named "cs" and "PDBOUT" within this folder, create them before 
                running Glutton script.

        Step3. Put the input chemical shift file (in NMRSTAR 3.0 format) in the subfolder "cs"

        Step4. change the input.txt to set parameters to determine the candidate pool

                bmr15398.str   # name of the input chemical shift file 
                2              # LEVEL of the database, 1 - high-resolution; 2 - medium-resolution; 3 - low-resolution (recommended 2)
                200            # Number of structures to be generated (recommended 2500 or more)
                0.5            # the selected width of the chemical shift distribution used to derive statistical distributions (recommended 0.5)
                0              # (recommended 0) don't use log of the density distributions of Phi and Psi 
                               # angles as probability function to generate angles
                               # if set to 1, log(p) will be used to generate Phi and Psi angles

        NOTE: if there are insufficient data to generate distributions, the choice can be 
        [1]. use the low-resolution database - LEVEL3;            
        [2]. increase the selected width of the chemical shift distribution. 
 
        Step5. use the following command to run script:

                python Glutton_R1.0.py
 
Part 3 of 4. Python libraries needed to run Glutton script

        The Glutton script has been tested using Anaconda python distribution version 3.6. 
        The Anaconda python package can be downloaded for free at: https://www.anaconda.com/download
        
        The compatibilty of Glutton script with Python 2.7 has NOT been evaluted. If you are using the default python  
        in linux, you are most likely using python 2.7. You will have to manually sepcify the python 
        interpretor as:
        
                /usr/bin/python3.5 Glutton_R1.0.py
        
        If you use the default python in Linux and don't have pip command, you will need to install pip as the following:
        
                sudo apt install python3-pip (Ubuntu)
        
                sudo yum install python3-pip (centos/redhat) - not tested

        Python libraries needed to run Glutton script:

        [1]. Install Biopython

                conda install biopython (or use "pip install biopython")
  
        Biopython webpage: http://biopython.org/
  
        [2]. Install PeptideBuilder

                pip install PeptideBuilder

        PeptideBuilder webpage: https://github.com/mtien/PeptideBuilder
  
        Reference:
        M. Z. Tien, D. K. Sydykova, A. G. Meyer, C. O. Wilke (2013). PeptideBuilder:
        A simple Python library to generate model peptides. PeerJ 1:e80.

        [3]. Install nmrstarlib
        
                pip install nmrstarlib

        nmrstarlib webpage: https://github.com/MoseleyBioinformaticsLab/nmrstarlib        

        [4]. Install other python libraries such as numpy, pandas if you don't have them.
        
                pip install numpy pandas

Part 4 of 4. Protein Chemical Shift - Structure Database (Glutton database)

        All database files 

                LEVEL1_ALL_20180116.dat  (high-resolution)

                LEVEL2_ALL_20180116.dat  (medium-resolution)

                LEVEL3_ALL_20180116.dat  (low-resolution)

        are in csv format. 

        Each row represents a single residue, and the meaning of each column is as follows:

        column 1: residue type (one letter code)

        column 2: the chemical shift value of the H on the backbone N

        column 3: the chemical shift value of the backbone N

        column 4: the chemical shift value of the H on the backbone C_Alpha

        column 5: the chemical shift value of the backbone C_Alpha

        column 6: the chemical shift value of the backbone C'

        column 7: Secondary structure type of this residue

        column 8: Phi angle of this residue

        column 9: Psi angle of this residue

        NOTE: if the chemical shift value of an atom is not available, Glutton sets it as 9999.0. 
