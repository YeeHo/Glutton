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
import matplotlib.pyplot as plt
import sys

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

def wavg(avg, weight):
    """ http://stackoverflow.com/questions/10951341/pandas-dataframe-aggregate-function-using-multiple-columns
    In rare instance, we may not have weights, so just return the mean. Customize this if your business case
    should return otherwise.
    """
    if weight.sum()==0:
        return math.nan
    try:
        return (avg * weight).sum() / weight.sum()
    except ZeroDivisionError:
        return math.nan

def GetCandidate(inputCS, ResNum, df, cutoff, weight):
    """Differences between exp and database records."""
    ATOMS = ['N1','CA1','H1','HA1','Cp1']
    iMark = [ 0,   0,    0,    0,    0]

    df_New = df[df['RES1'] == inputCS['RES1'].iloc[ResNum]].copy()
    df_New[['H1', 'N1', 'HA1', 'CA1', 'Cp1']] = df_New[['H', 'N', \
       'HA', 'CA', 'Cp']] - inputCS[['H', 'N', 'HA', 'CA', 'Cp']].iloc[ResNum].values

    if abs(inputCS['N'].iloc[ResNum]) < 500.0:
        cutoff[0] = df_New[(df_New['N1']<500.0)&(df_New['N']<500.0)]['N'].quantile(0.90) - df_New[(df_New['N1']<500.0)&(df_New['N']<500.0)]['N'].quantile(0.1)
        df_New = df_New[(df_New['N1'].abs() < (cutoff[0]*weight))&(df_New['N']<500.0)]
        iMark[0]=1
    if abs(inputCS['CA'].iloc[ResNum]) < 500.0:
        cutoff[1] = df_New[(df_New['CA1']<500.0)&(df_New['CA']<500.0)]['CA'].quantile(0.90) - df_New[(df_New['CA1']<500.0)&(df_New['CA']<500.0)]['CA'].quantile(0.1)
        df_New = df_New[(df_New['CA1'].abs() < (cutoff[1]*weight))&(df_New['CA']<500.0)]
        iMark[1]=1
    if abs(inputCS['H'].iloc[ResNum]) < 500.0:
        cutoff[2] = df_New[(df_New['H1']<500.0)&(df_New['H']<500.0)]['H'].quantile(0.90) - df_New[(df_New['H1']<500.0)&(df_New['H']<500.0)]['H'].quantile(0.1)
        df_New = df_New[(df_New['H1'].abs() < (cutoff[2]*weight))&(df_New['H']<500.0)]
        iMark[2]=1
    elif abs(inputCS['HA'].iloc[ResNum]) < 500.0:
        cutoff[3] = df_New[(df_New['HA1']<500.0)&(df_New['HA']<500.0)]['HA'].quantile(0.90) - df_New[(df_New['HA1']<500.0)&(df_New['HA']<500.0)]['HA'].quantile(0.1)
        df_New = df_New[(df_New['HA1'].abs() < (cutoff[3]*weight))&(df_New['HA']<500.0)]
        iMark[3]=1
    elif abs(inputCS['Cp'].iloc[ResNum]) < 500.0:
        cutoff[4] = df_New[(df_New['Cp1']<500.0)&(df_New['Cp']<500.0)]['Cp'].quantile(0.90) - df_New[(df_New['Cp1']<500.0)&(df_New['Cp']<500.0)]['Cp'].quantile(0.1)
        df_New = df_New[(df_New['Cp1'].abs() < (cutoff[4]*weight))&(df_New['Cp']<500.0)]
        iMark[4]=1

    AvaAtom = []
    DiffAvaAtom = []
    iSign = []
    cutoffX = []
    for i in range(5):
        if iMark[i] == 1:
            AvaAtom.append(ATOMS[i])
            xAtomS = ATOMS[i]+'S'
            df_New[xAtomS] = df_New[ATOMS[i]] - df_New[ATOMS[i]].sum()
            DiffAvaAtom.append(xAtomS)
            if df_New[ATOMS[i]].sum() > 0:
                iSign.append(False)
            else:
                iSign.append(True)
            cutoffX.append(cutoff[i]*weight)

    for i in range(int(0.8*len(df_New.index))):
        for x in range(len(AvaAtom)):
            df_New[DiffAvaAtom[x]] = df_New[AvaAtom[x]] - df_New[AvaAtom[x]].sum()
            if df_New[AvaAtom[x]].sum() > 0:
                iSign[x] = False
            else:
                iSign[x] = True

        df_TMP = df_New[DiffAvaAtom].copy()
        df_New['Tot'] = df_TMP[DiffAvaAtom].pow(2).sum(axis=1)

        df_New = df_New.sort_values('Tot', ascending=False)
        for x in range(len(AvaAtom)):
            for y in range(x+1, len(AvaAtom)):
                if df_New[DiffAvaAtom[x]].abs().iloc[0] < df_New[DiffAvaAtom[y]].abs().iloc[0]:
                    tmp = DiffAvaAtom[x]
                    DiffAvaAtom[x] = DiffAvaAtom[y]
                    DiffAvaAtom[y] = tmp
                    boolX = iSign[x]
                    iSign[x] = iSign[y]
                    iSign[y] = boolX

        df_New = df_New.sort_values(['Tot']+DiffAvaAtom, ascending=[True]+iSign)
        df_New.drop(df_New.index[0], inplace=True)

        if len(df_New.index) < 72:
            break

    return df_New, len(df_New.index), inputCS['RES1'].iloc[ResNum]

def GetPhiPsiCandidate(seq, inputPhiPsi, ResNum, df, cutoff):
    """Differences between exp and database records."""

    df_New = df[df['RES1'] == seq[ResNum]]
    # print(inputPhiPsi[['PHI','PSI']].iloc[ResNum].values)
    with warnings.catch_warnings():
         warnings.simplefilter("ignore")
         df_New[['PHI1','PSI1']] = df_New[['PHI','PSI']] - inputPhiPsi[['PHI','PSI']].iloc[ResNum].values
         df_New[['PHI1','PSI1']] = df_New[['PHI1','PSI1']].abs()

    df_New = df_New[df_New['PHI1'] < cutoff]
    df_New = df_New[df_New['PSI1'] < cutoff]

    x = plt.hist(df_New[df_New['H']<100.0 ]['H'],bins=50)
    df_H = wavg(x[1][0:50],x[0])
    x = plt.hist(df_New[df_New['N']<250.0]['N'],bins=50)
    df_N = wavg(x[1][0:50],x[0])
    x = plt.hist(df_New[df_New['HA']<100.0]['HA'],bins=50)
    df_HA = wavg(x[1][0:50],x[0])
    x = plt.hist(df_New[df_New['CA']<250.0]['CA'],bins=50)
    df_CA = wavg(x[1][0:50],x[0])
    x = plt.hist(df_New[df_New['Cp']<250.0]['Cp'],bins=50)
    df_Cp = wavg(x[1][0:50],x[0])

    return df_H, df_N, df_HA, df_CA, df_Cp

def AnglesCDF(data, NumofStructGen):
    """"Calculate CDF for data."""
    data = data[data<359.0]
    hist, bins = np.histogram(data, bins=72)
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
def build_phi_psi_model(pdb_filename):
    parser=PDBParser()
    structure=parser.get_structure('sample', \
                                    path.join('./PDB/', pdb_filename))
    model=structure[0]
    chain=model['A']
    seq=""
    phi_diangle=[]
    psi_diangle=[]
    omega_diangle=[]
    for res in chain:
        if(res.get_resname() in resdict.keys()):

            seq+=resdict[res.get_resname()]
            if(len(seq)==1):
                N_prev=res['N']
                CA_prev=res['CA']
                C_prev=res['C']
            else:
                n1=N_prev.get_vector()
                ca1=CA_prev.get_vector()
                c1=C_prev.get_vector()

                C_curr=res['C']
                N_curr=res['N']
                CA_curr=res['CA']

                c=C_curr.get_vector()
                n=N_curr.get_vector()
                ca=CA_curr.get_vector()

                psi= calc_dihedral(n1, ca1, c1, n) ##goes to current res
                omega= calc_dihedral(ca1, c1, n, ca)
                phi= calc_dihedral(c1, n, ca, c) ##goes to current res

                phi_diangle.append(phi*180.0/math.pi)
                psi_diangle.append(psi*180.0/math.pi)
                # omega_diangle.append(omega*180.0/math.pi)

                N_prev=res['N']
                CA_prev=res['CA']
                C_prev=res['C']

    return seq, phi_diangle, psi_diangle

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

        dfxPhi = dfx.copy()
        dfxPsi = dfx.copy()

        Phi[inum,:] = AnglesCDF(dfxPhi['PHI'].values, NumofStructGen)
        Psi[inum,:] = AnglesCDF(dfxPsi['PSI'].values, NumofStructGen)

    return seq, Phi, Psi, NumofStructGen

def GenerateCS(PDBFile, cutoff, databaseLevel):
    """"Generate Phi and Psi angle ensembles based on statistics."""
    seq, phi_d, psi_d = build_phi_psi_model(PDBFile)

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
    seqLen = len(seq)
    phi_diangle = np.asarray([360.0] + phi_d)
    psi_diangle = np.asarray(psi_d + [360.0])

    PhiPsi = pd.DataFrame(data=np.column_stack((phi_diangle, psi_diangle)), \
         columns=['PHI','PSI'], dtype=np.float)

    CSOUT = np.zeros(seqLen*5).reshape(5,seqLen)
    seqList = list(seq)
    for inum in range(seqLen):
        df_H, df_N, df_HA, df_CA, df_Cp = GetPhiPsiCandidate(seqList, PhiPsi, inum, dfSS, cutoff)
        CSOUT[0,inum] = df_H
        CSOUT[1,inum] = df_N
        CSOUT[2,inum] = df_HA
        CSOUT[3,inum] = df_CA
        CSOUT[4,inum] = df_Cp

    return CSOUT

def constructStr(seq, Phi, Psi, NumofStructGen):
    """"Construct structures based on given Phi and Psi angle ensembles."""
    SeqLen = len(Phi[:,0])
    iCount = 0;
    PhiOut = np.zeros((int(SeqLen), int(NumofStructGen)))
    PsiOut = np.zeros((int(SeqLen), int(NumofStructGen)))
    for j in range(10*NumofStructGen):
        model_structure_phi_psi= PeptideBuilder.make_structure(seq, Phi[:,j], Psi[:,j])
        Clashes = disCheck(model_structure_phi_psi)
        if Clashes == False:
            PhiOut[:,iCount] = Phi[:,j]
            PsiOut[:,iCount] = Psi[:,j]

            iCount = iCount + 1
            if iCount % 100 == 0:
               print(iCount)
            out = PDBIO()
            out.set_structure(model_structure_phi_psi)
            name = 'output'+str(iCount)+'.pdb'
            out.save('./PDBOUT/'+name)
        if iCount > (NumofStructGen-1):
            break
    return PhiOut, PsiOut

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
    for i in range(6):
        boards.append(f.readline().strip())
    return boards

def glutton_cs2str(inputPara):
    """"Generate structure ensemble based on the input.txt."""

    print("Number of the generated structure will be updated every 100")
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

        PhiOutX, PsiOutX = constructStr(seqncbd, NCBDPhi, NCBDPsi, int(inputPara[2]))
    else:
        PhiOutX, PsiOutX = constructStr(seq, Phi, Psi, int(inputPara[2]))

    PhiOut = np.asarray(Phi)
    PsiOut = np.asarray(Psi)

    if int(inputPara[4]) == 0:
        np.savetxt('PhiAngles.txt', PhiOut, fmt='%7.2f')
        np.savetxt('PsiAngles.txt', PsiOut, fmt='%7.2f')
    elif int(inputPara[4]) == 1:
        np.savetxt(inputPara[0]+'PhiAngles.txt', PhiOutX, fmt='%7.2f')
        np.savetxt(inputPara[0]+'PsiAngles.txt', PsiOutX, fmt='%7.2f')
        np.savetxt(inputPara[0]+'PhiAnglesTot.txt', PhiOut, fmt='%7.2f')
        np.savetxt(inputPara[0]+'PsiAnglesTot.txt', PsiOut, fmt='%7.2f')

def glutton_str2cs(inputPara):
    """"Generate chemical shift distributions based on the input.txt."""

    CSOUT = GenerateCS(inputPara[0], float(inputPara[3]), int(inputPara[1]))
    seq, phi_d, psi_d = build_phi_psi_model(inputPara[0])
    SeqLen = len(seq)
    seqList = list(seq)

    print('{0:>6s} {1:>3s} {2:>6s} {3:>6s} {4:>6s} {5:>6s} {6:>6s}'.format( \
        'RESNUM', 'RES', 'H', 'N', 'HA', 'CA', "C'"))
    for i in range(SeqLen):
        print('{0:>3d} {1:>2s} {2:>6.2f} {3:>6.2f} {4:>6.2f} {5:>6.2f} {6:>6.2f} \
              '.format(i+1, seqList[i], CSOUT[0,i], CSOUT[1,i], CSOUT[2,i], CSOUT[3,i], CSOUT[4,i],))


def glutton():
    inputPara = readinput()
    print("Input Parameters: ", inputPara)

    if int(inputPara[5]) == 0:
        glutton_str2cs(inputPara)
    else:
        glutton_cs2str(inputPara)

glutton()
