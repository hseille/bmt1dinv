'''
    File name: tree.py
    Authors: Hoël Seillé / Gerhard Visser
    Date created: 01/06/2020
    Date last modified: 26/10/2020
    Python Version: 3.6
'''
__author__ = "Hoël Seillé / Gerhard Visser"
__copyright__ = "Copyright 2020, CSIRO"
__credits__ = ["Hoël Seillé / Gerhard Visser"]
__license__ = "GPLv3"
__version__ = "1.0.0"
__maintainer__ = "Hoël Seillé"
__email__ = "hoel.seille@csiro.au"
__status__ = "Beta"



import pandas as pd
import numpy as np
from scipy import linalg


class Dim_Tree():
    def __init__(self, con, left=None, right=None, p=None, isLeaf=False):
        self.con = con   # condition / attribute upon the tree splits (None when Leaf)
        self.left  = left   # left is for True condition (None when Leaf)
        self.right = right   # right is for False condition (None when Leaf)
        self.p = p   # (None when Split)
        self.isLeaf = isLeaf   # boolean for Leaf of Split node
        
    def vg(self, exp, rho):  # variogram parameters
        self.vgExp = exp
        self.vgRho = rho
        
    def atts(self,atts_list):  # attributes list
        self.atts = atts_list



def readTree(tree):
    global nr
    if tree['node'][nr] == 'Leaf':
        p = tree['val'][nr]
        return Dim_Tree(None, None,None,p,True)
    else:
        con = tree['val'][nr]
        nr += 1
        left = readTree(tree)
        nr += 1
        right = readTree(tree)
        return Dim_Tree(con, left, right , None, False)



def readTreeFile(tree_file, atts_dic):
    f=open(tree_file, "r")
    f1=f.readlines()
    # store the variogram parameters
    vgExp = float(f1[0])
    vgRho = float(f1[1])
    
    #read the .txt file
    tree_txt = pd.read_csv(tree_file,header=None,skiprows=2, 
                       delim_whitespace=True, 
                       names=['depth','node','val'])
    global nr
    nr = 0
    #read and store the decision tree
    tr = readTree(tree_txt)
    tr.vg(vgExp, vgRho)
    tr.atts(list(atts_dic))
    
    return tr 


def readAttsFile(atts_file):
    #read the .txt file
    with open(atts_file, "r") as f:
        atts_dic=[]
        for line in f:
            atts_dic.append(line.split(','))

    return atts_dic


def attMap(data,atts_dic):
    #Define attributes to use
    elip = data['elipM']
    beta = data['betaM']
    difPol = data['difPolM']
    attDef = []
    for i in range(len(atts_dic)):
        attDef.append([eval(atts_dic[i][0]), 
                       eval(atts_dic[i][1]), 
                       eval(atts_dic[i][2])])
    atts = [attDef[i][2] for i in range(len(attDef))]
    return atts



def leaffor(tr, atts):
    if tr.isLeaf:
        return tr.p
    else:
        con = int(tr.con)
        if atts[con-1]:
            tr = tr.left
            return leaffor(tr, atts)
        else:
            tr = tr.right
            return leaffor(tr, atts)



def covar(d, vgExp, vgRho):
    v = d/vgRho
    v = v**vgExp
    return np.exp(-v)



def getC(data, tr, nsCov=None, noiseCovar=False):

    # noiseCovar: in case a full covariance matrix is available from the 
    # processing error. If True, nsCov needs to be provided
	
    nF = len(data)
    attss = attMap(data, tr.atts)
    attss = np.array(attss).T
    ps = [leaffor(tr, attss[i,:]) for i in range(nF)]
    ps = np.array(ps)

    if noiseCovar:
        ns = np.diag(nsCov)**0.5
        ss = (ps**2 + ns **2)**0.5
    else:
        ns = data['ZdetLnSd']  
        ss = (ps**2 + ns **2)**0.5
        
    C = np.zeros((nF,nF))
    dlf = np.zeros((nF,nF)) # log distance between frequencies for the variogram
    for i in range(nF):
        for j in range(nF):
            dlf[i,j] = abs(np.log10(data['freq'][i])-np.log10(data['freq'][j]))
            C[i,j] = covar(dlf[i,j],tr.vgExp, tr.vgRho)
            C[i,j] *= ps[i]*ps[j]
            
    if noiseCovar:
        C += ns**2
    else:
        for i in range(nF):
            C[i,i] += ns[i]**2
    return C, ss



def exportForTransD(siteId,data, tr, errorfloor=-1, fcorr = False, min_errorfloor = 0.01):
    nF = len(data)
    if errorfloor > 0:
        C = np.zeros((nF,nF))
        ss = np.zeros((nF))
        for i in range(nF):
            if np.log(1+errorfloor) > data['ZdetLnSd'][i]:
                C[i,i] = (np.log(1+errorfloor))**2
            else:
                C[i,i] = (data['ZdetLnSd'][i])**2
            ss[i] = pd.Series(C[i,i]**0.5)
        
    else:
        C, ss = getC(data, tr)
        # we use a minimum error floor of 1% in case the overall error is too small
        # min_errorfloor = 0.01
        for i in range(nF):
            if np.log(1+min_errorfloor) > ss[i]:
                C[i,i] = (np.log(1+min_errorfloor))**2
            ss[i] = pd.Series(C[i,i]**0.5)
			
			

    if fcorr:
        # Calculate inverse of C
        # D = linalg.inv(C)
        D = np.linalg.lstsq(C, np.identity(C.shape[0]))
        print(D.shape)
  
    else:
        D = np.diag(ss**-2)
		
    #write to csv file   
    ss=pd.Series(ss)
    df = pd.concat([data['freq'], data['ZdetRLn'], data['ZdetILn'], ss], axis=1)
    df = pd.concat([df, pd.DataFrame(D)], axis=1)
    columns = ['freq','Zr','Zi','std'] + ['D%i'%(i) for i in range(1,nF+1)]
    df.columns = columns
	
	return df,ss


