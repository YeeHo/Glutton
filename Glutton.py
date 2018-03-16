import pandas as pd
import numpy as np
import nmrstarlib
from Bio.PDB import PDBParser
from Bio.PDB import Superimposer
from Bio.PDB.Atom import *
from Bio.PDB.Residue import *
from Bio.PDB.Chain import *
from Bio.PDB.Model import *
from Bio.PDB.Structure import *
from Bio.PDB.Vector import *
from Bio.PDB.Entity import*
import math
from PeptideBuilder import Geometry
import PeptideBuilder
from os import path
from Bio.PDB import *

resdict = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
           'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
           'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
           'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}

PDBdir = "./pdb/"

##########################################################################
# Section 1: Get Chemical shift data from BMRB File (NMRSTAR 3.1 format) #
##########################################################################

def LoadNMRCS(BMRBID,FOLDER):
    """Load a file in NMRSTAR 3.1 format and get Chemical shift data."""
    # print(FOLDER+BMRBID)
    starfiles = nmrstarlib.read_files(FOLDER+BMRBID)
    starfile = starfile = list(starfiles)
    keyList = list(starfile[0].keys())
    indices = [i for i, s in enumerate(keyList) if 'comment' in s]
    df, iNum = getCS(starfile, keyList, indices)
    return df, iNum

def LoadPDB(PDBID,FOLDER):
    from Bio.PDB import PDBParser
    pdbfile = FOLDER+PDBID+'.pdb'
    parser = PDBParser()
    structure = parser.get_structure(PDBID, pdbfile)
    return structure

def getCS(starfile, keyList, indices):
    """read Chemical shift data."""
    ThreeLetterCode = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY',
                       'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                       'THR', 'TRP', 'TYR', 'VAL']

    OneLetterCode = ['A',   'R',   'N',   'D',   'C',   'E',   'Q',   'G',
                     'H',   'I',   'L',   'K',   'M',   'F',   'P',   'S',
                     'T',   'W',   'Y',   'V']

    ATOM_List = ['N', 'H', 'HA', 'HA2', 'HA3', 'HB', 'HB2',
                 'HB3', 'CA', 'CB', 'C']

    iResNum = 'XXX'
    iNum = 0
    df_List = [['999.0' for x in range(12)] for y in range(1000)]
    for i in range(0,len(indices)):
        if ('Chemical Shift Ambiguity Index Value Definitions' in starfile[0][keyList[indices[i]]]):
            fname_cs = keyList[indices[i]+1]
            loop_name = getLoopCS(starfile,fname_cs)
            if loop_name == 'loop_x' or loop_name == 'loop_49':
                print('Please Check the Format of Your Chemical Shift File!')
            else:
                ind_len = len(starfile[0][fname_cs][loop_name][1])
                ShortestID = starfile[0]['data']

                for ind in range(1, ind_len):
                    bmrdid = starfile[0]['data']
                    tComp_index_ID = starfile[0][fname_cs][loop_name][1][ind-1]["Atom_chem_shift.Comp_index_ID"]
                    tCS_Val = starfile[0][fname_cs][loop_name][1][ind-1]["Atom_chem_shift.Val"]
                    tAsym_ID = starfile[0][fname_cs][loop_name][1][ind-1]["Atom_chem_shift.Auth_seq_ID"]
                    tSeq_ID = starfile[0][fname_cs][loop_name][1][ind-1]["Atom_chem_shift.Seq_ID"]
                    tComp_ID = starfile[0][fname_cs][loop_name][1][ind-1]["Atom_chem_shift.Comp_ID"]
                    tAtom_ID = starfile[0][fname_cs][loop_name][1][ind-1]["Atom_chem_shift.Atom_ID"]
                    tAssem_ID = starfile[0][fname_cs][loop_name][1][ind-1]["Atom_chem_shift.Entity_assembly_ID"]
                    tEntity_ID = starfile[0][fname_cs][loop_name][1][ind-1]["Atom_chem_shift.Entity_ID"]

                    tOneLetter = tComp_ID
                    for i in range(0,20):
                        tOneLetter = tOneLetter.replace(ThreeLetterCode[i], OneLetterCode[i])

                    if tComp_ID.strip() in ThreeLetterCode:
                       if iResNum == 'XXX':
                          for xx in range(5,12):
                              df_List[iNum][xx] == '999.0'
                          df_List[iNum][0] = bmrdid
                          df_List[iNum][1] = tComp_index_ID
                          df_List[iNum][2] = tAsym_ID
                          df_List[iNum][3] = tComp_ID
                          df_List[iNum][4] = tOneLetter
                          df_List = addCStoFrame(df_List,iNum,tAtom_ID, tCS_Val)
                          iResNum = tComp_index_ID
                       else:
                          if iResNum != tComp_index_ID:
                              iNum = iNum + 1
                              for xx in range(5,12):
                                  df_List[iNum][xx] == '999.0'
                              df_List[iNum][0] = bmrdid
                              df_List[iNum][1] = tComp_index_ID
                              df_List[iNum][2] = tAsym_ID
                              df_List[iNum][3] = tComp_ID
                              df_List[iNum][4] = tOneLetter
                              df_List = addCStoFrame(df_List,iNum,tAtom_ID, tCS_Val)
                              iResNum = tComp_index_ID
                          else:
                              df_List[iNum][0] = bmrdid
                              df_List[iNum][1] = tComp_index_ID
                              df_List[iNum][2] = tAsym_ID
                              df_List[iNum][3] = tComp_ID
                              df_List[iNum][4] = tOneLetter
                              df_List = addCStoFrame(df_List,iNum,tAtom_ID, tCS_Val)

    df = pd.DataFrame(df_List[0:iNum+1][:],columns=['BMRBID','RESNUMCS','RESNUMPDB','RES3','RES1',
                                                    'H','N','HA','HB','CA','CB','Cp'])
    return df, iNum

def getLoopCS(starfile,fname_cs):
    loop_name = 'loop_x'
    for i in range(50):
        loopTest = 'loop_'+str(i)
        try:
            if len(starfile[0][fname_cs][loopTest][0])>10:
               return loopTest
        except:
            pass
    return loop_name

def addCStoFrame(df,iNum,tAtom_ID, tCS_Val):
    if tAtom_ID == 'H':
        df[iNum][5] = tCS_Val
    elif tAtom_ID == 'N':
        df[iNum][6] = tCS_Val
    elif tAtom_ID == 'HA' or tAtom_ID == 'HA2' or tAtom_ID == 'HA3':
        df[iNum][7] = tCS_Val
    elif tAtom_ID == 'HB' or tAtom_ID == 'HB3' or tAtom_ID == 'HB3':
        df[iNum][8] = tCS_Val
    elif tAtom_ID == 'CA':
        df[iNum][9] = tCS_Val
    elif tAtom_ID == 'CB':
        df[iNum][10] = tCS_Val
    elif tAtom_ID == 'C':
        df[iNum][11] = tCS_Val
    return df

#######################################################
# Section 2: Load Database and identify distributions #
#######################################################

def LoadData(FileName):
    """Load Database."""
    names = ['RES1','H','N','HA','CA','Cp','SS','PHI', 'PSI']
    dtype = {'RES1':object, 'H':float, 'N':float, 'HA':float, 'CA':float,
             'Cp':float, 'SS':object, 'PHI':float, 'PSI':float}
    return pd.read_csv(FileName, names=names, dtype=dtype)

def GetSTD(df, ResID):
    """Calculate STD for a Dataframe."""
    ATOMS = ['H','N','HA','CA','Cp']
    dfstd = [0.0, 0.0, 0.0, 0.0, 0.0]
    df = df[df['RES1'] == ResID]
    for i in range(len(ATOMS)):
        dx = df[df[ATOMS[i]]<500.0]
        dfstd[i] = dx[ATOMS[i]].std()
    return dfstd

def GetCandidate(inputCS, ResNum, df, cutoff, weight):
    """Differences between exp and database records."""
    ATOMS = ['H1','N1','HA1','CA1','Cp1']
    df_New = df[df['RES1'] == inputCS['RES1'].iloc[ResNum]].copy()
    df_New[['H1', 'N1', 'HA1', 'CA1', 'Cp1']] = df_New[['H', 'N', \
       'HA', 'CA', 'Cp']] - inputCS[['H', 'N', 'HA', 'CA', 'Cp']].iloc[ResNum].values

    df_New[['H1', 'N1', 'HA1','CA1', 'Cp1']] = df_New[['H1', 'N1', 'HA1', 'CA1', 'Cp1']].abs()

    if abs(inputCS['N'].iloc[ResNum]) < 500.0:
        df_New = df_New[df_New['N1'] < (cutoff[1]*weight)]
    if abs(inputCS['CA'].iloc[ResNum]) < 500.0:
        df_New = df_New[df_New['CA1'] < (cutoff[3]*weight)]
    if abs(inputCS['H'].iloc[ResNum]) < 500.0:
        df_New = df_New[df_New['H1'] < (cutoff[0]*weight)]
    elif abs(inputCS['HA'].iloc[ResNum]) < 500.0:
        df_New = df_New[df_New['HA1'] < (cutoff[2]*weight)]
    elif abs(inputCS['Cp'].iloc[ResNum]) < 500.0:
        df_New = df_New[df_New['Cp1'] < (cutoff[4]*weight)]

    return df_New, len(df_New.index), inputCS['RES1'].iloc[ResNum]

def AnglesCDF(data, NumofStructGen):
    """"Calculate CDF for data."""
    hist, bins = np.histogram(data, bins=36)
    hist[hist==0] = 1
    bin_midpoints = bins[:-1] + np.diff(bins)/2
    cdf = np.cumsum(hist)
#    print(cdf[-1])
    cdf = cdf / cdf[-1]
    values = np.random.rand(10*NumofStructGen)
    value_bins = np.searchsorted(cdf, values)
    angles = bin_midpoints[value_bins]
    return angles

#######################################################
#  Section 3: Generate structures from distributions  #
#######################################################

def GeneratePhiPsi(CSFile, databaseLevel, NumofStructGen, Weight):
    """"Generate Phi and Psi angle ensembles based on statistics."""
    dx, iLen = LoadNMRCS(CSFile, './cs/')
    with warnings.catch_warnings():
         warnings.simplefilter("ignore")
         inputCSdf = dx.convert_objects(convert_numeric=True)
    NumList = dx['RESNUMCS'].values
    seq = dx['RES1'].values

    if databaseLevel == 1:
        df = LoadData('LEVEL1_ALL_20180116.dat')
    elif databaseLevel == 2:
        df = LoadData('LEVEL2_ALL_20180116.dat')
    elif databaseLevel == 3:
        df = LoadData('LEVEL3_ALL_20180116.dat')
    else:
        print('Please input a number among 1, 2, 3 to select among  \
                three precision levels:')
        print('1. High resolution database')
        print('2. Medium resolution database')
        print('3. Low resolution database')

    dfSS = df.copy()

    Phi = np.random.rand(10*iLen*NumofStructGen).reshape(iLen,10*NumofStructGen,)
    Psi = np.random.rand(10*iLen*NumofStructGen).reshape(iLen,10*NumofStructGen)

    for inum in range(iLen):
        cutoff = GetSTD(df,inputCSdf['RES1'].iloc[inum])
        dfx, xLen, ResID = GetCandidate(inputCSdf, inum, dfSS, cutoff, Weight)

        Phi[inum,:] = AnglesCDF(dfx['PHI'].values, NumofStructGen)
        Psi[inum,:] = AnglesCDF(dfx['PSI'].values, NumofStructGen)

    return seq, Phi, Psi, NumofStructGen

def constructStr(seq, Phi, Psi, NumofStructGen):
    """"Construct structures based on given Phi and Psi angle ensembles."""
    iCount = 0;
    for j in range(10*NumofStructGen):
        model_structure_phi_psi= PeptideBuilder.make_structure(seq, Phi[:,j], Psi[:,j])
        Clashes = disCheck(model_structure_phi_psi)
        if Clashes == False:
            iCount = iCount + 1
            if iCount % 100 == 0:
               print(iCount)
            out = PDBIO()
            out.set_structure(model_structure_phi_psi)
            name = 'output'+str(iCount)+'.pdb'
            out.save('./PDBOUT/'+name)
        if iCount > (NumofStructGen-1):
            break

def disCheck(structure):
    """"Check Ca overlap."""
    model = structure[0]
    chain = model['A']
    for residue1 in chain:
        for residue2 in chain:
            if residue1 != residue2:
                distance = residue1['CA'] - residue2['CA']
                if distance < 1.5:
                    return True
    return False

########################################################################
#         Section 4: Construct structures for selected proteins        #
########################################################################

def readinput():
    """"Read input parameters."""
    f = open('input.txt')
    boards = []
    for i in range(5):
        boards.append(f.readline().strip())
    return boards

def glutton():
    """"Generate structure ensemble based on the input.txt."""
    inputPara = readinput()
    print("Input Parameters: ", inputPara)
    print("Number of the generated tructure will be updated every 100")
    seq, Phi, Psi, NumofStructGen = GeneratePhiPsi(inputPara[0], int(inputPara[1]),
            int(inputPara[2]), float(inputPara[3]))

    if inputPara[0] == 'bmr15398.str':
        print('NCBD test activated!')
        NCBDPhi = np.random.rand(540*int(inputPara[2])).reshape(54,10*int(inputPara[2]))
        NCBDPsi = np.random.rand(540*int(inputPara[2])).reshape(54,10*int(inputPara[2]))
        seqncbd = 'ISPSALQDLLRTLKSPSSPQQQQQVLNILKSNPQLMAAFIKQRTAKYVANQPGMQ'

        # Because of the missing chemical shift data for two residues
        NCBDPhi[0:17, :] = Phi[0:17,:]
        NCBDPhi[17:18,:] = [-62.2154587855]*10*int(inputPara[2])
        NCBDPhi[18:42,:] = Phi[17:41,:]
        NCBDPhi[42:43,:] = [66.4422232256]*10*int(inputPara[2])
        NCBDPhi[43:,:] = Phi[41:,:]

        NCBDPsi[0:17, :] = Psi[0:17,:]
        NCBDPsi[17:18,:] = [-21.0976186348]*10*int(inputPara[2])
        NCBDPsi[18:42,:] = Psi[17:41,:]
        NCBDPsi[42:43,:] = [14.3591021712]*10*int(inputPara[2])
        NCBDPsi[43:,:] = Psi[41:,:]

        constructStr(seqncbd, NCBDPhi, NCBDPsi, int(inputPara[2]))
    else:
        constructStr(seq, Phi, Psi, int(inputPara[2]))

glutton()
