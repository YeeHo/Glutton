![alt text](https://github.com/YeeHo/Glutton/blob/master/NCBD-49.png)
Figure 1. Side-by-side comparison of the Phi and Psi angle distributions based on 10,000 structures obtained from 10-us MD simulations(left) and Glutton (right) for ALA-49 in NCBD. 

# Glutton
Please take 5 minutes to read this brief readme file before running the Glutton script. 

Part 1 of 4. Introduction of Glutton

        The Glutton database and a complementary python script are developed for constructing structure 
        ensembles solely based on backbone chemical shift data for flexible proteins.
        
        The Glutton database that we have developed is specifically designed to fill a void
        currently existing in the bioinformatics tools for the structural analysis of intrinsically
        disordered proteins (IDPs) based on experimental information, namely Nuclear Magnetic Resonance (NMR). 
        Our Glutton database contains both chemical shifts and their corresponding structural information 
        for a total of 5,270 proteins. The Glutton database is organized hierarchically in three tiers
        (levels of resolution) that permit customization of the database for specific applications: 
        for calculations of structural ensembles of IDPs, or for determination and prediction of well-folded 
        protein structures. And, most importantly, Glutton interprets chemical shifts as ensemble 
        averages, and thus calculates conformational distributions from chemical shift datasets. 
        These distributions permit the straightforward determination of a maximally broad structural 
        ensemble consistent with the data, and thus it optimizes sampling of the conformational space that is
        available to the IDP as shown in the above plot. These ensembles provide a realistic description
        of the conformational heterogeneity and structural propensities of IDPs, but they can also be used 
        as starting points for the refinement of protein structure determination (adding additional 
        experimental information) or protein prediction (incorporating a molecular force-field).

Part 2 of 4. Use Glutton script to generate structure ensembles based on chemical shifts or vice versa.

        Step1. Put database files (LEVEL*.dat), Glutton.py and input.txt files in the same folder.

        Step2. If there are no subfolders named "cs" and "PDBOUT" within this folder, create them before 
                running the Glutton script.

        Step3. Put the input chemical shift file (in NMRSTAR 3.0 format) in the subfolder "cs" or pdb file in subfolder "pdb"

        Step4. Change the input.txt to set parameters to determine the characteristics of the structure 
        ensemble to be generated.

                bmr15398.str   # name of the input chemical shift file in cs folder
                2              # LEVEL of the database, 1 - high-resolution; 2 - medium-resolution; 3 - low-resolution 
                200            # Number of structures to be generated (when predicting cs from a structure, set it to 1)
                0.5            # the selected width of the chemical shift distribution to derive statistical distributions (or the phi/psi angle range if predicting cs from a structure)
                1              # 0 - Output all Phi and Psi angles before size exclusion; 1 - output Phi and Psi angles used in all the outputed structures
                1              # 0 - predict chemical shifts based on PDB structure; 1 - predict strcutures based on chemical shifts

        NOTE: if the available data are insufficient to generate distributions, you can either: 
        [1]. use the low-resolution database - LEVEL3;            
        [2]. increase the selected width of the chemical shift distribution. 
        [3]. Output Phi and Psi angle files are named "PhiAngles.txt" and "PsiAngles.txt". Each row contains all
             the generated angles corresponding to one residue.
 
        Step5. Use the following command to run the Glutton script:

                python Glutton.py

        NOTE: Example input files for Cs to structure is "input_cs2str_ncbd.txt" and for structure to CS is
         "input_str2cs_ncbd.txt". To use anyone of them as input file for Glutton, please change the file name
         to "input.txt" and put in the same folder as Glutton.py.

Part 3 of 4. Python libraries needed to run this script

        The Glutton script has been tested using the Anaconda python distribution version 3.6. 
        The Anaconda python package can be downloaded for free from: https://www.anaconda.com/download
        
        The compatibility of the Glutton script with Python version 2.7 has NOT been evaluated. If you are using the 
        default python in linux, you are most likely using python version 2.7. You will have to manually specify 
        the python interpreter as:
        
                /usr/bin/python3.5 Glutton.py
        
        If you use the default python in Linux and don't have the pip command, you will need to install pip 
        as follows:
        
                sudo apt install python3-pip (Ubuntu)
        
                sudo yum install python3-pip (centos/redhat) - not tested

        The python libraries needed to run Glutton script are:

        [1]. Install Biopython

                conda install biopython (or use "pip install biopython")
  
        Biopython webpage: http://biopython.org/
  
        [2]. Install PeptideBuilder

                pip3 install PeptideBuilder

        PeptideBuilder webpage: https://github.com/mtien/PeptideBuilder
  
        Reference:
        M. Z. Tien, D. K. Sydykova, A. G. Meyer, C. O. Wilke (2013). PeptideBuilder:
        A simple Python library to generate model peptides. PeerJ 1:e80.

        [3]. Install nmrstarlib
        
                pip3 install nmrstarlib

        nmrstarlib webpage: https://github.com/MoseleyBioinformaticsLab/nmrstarlib        

        [4]. Install other python libraries such as numpy, pandas if you don't have them.
        
                pip3 install numpy pandas

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
