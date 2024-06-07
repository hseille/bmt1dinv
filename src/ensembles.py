'''
    File name: ensembles.py
    Authors: Hoël Seillé / Gerhard Visser
    Date created: 01/10/2020
    Date last modified: 14/04/2021
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



import numpy as np
import struct
import pandas as pd
from numba import jit
import math
import random

import warnings
warnings.filterwarnings('ignore')

import MT as MT



"""
Read outputs of the MCMC trans-dimensional algorithm

Input:
Binary files
- ensemble of models
- ensemble of responses
- input csv file (observed data) Z in Log 

Output:
- ensembles
- input data Z in Linear scale (ohm)

"""

def sampModels(samp_models_file):

    f = open(samp_models_file, "rb")
    
    lkh = []
    prior =[]
    nLayers = []
    ensemble_models = []
    
    try:
        while True:
            byte = f.read(4)
            if not byte:
                break
            else:
                a = struct.unpack('=i', byte)[0]
                if a > 100000 or a < 0:
                    pass
                else:
                    nLayers.append(a)
                    byte = f.read(8)
                    lkh.append(struct.unpack('=d', byte)[0])
                    byte = f.read(8)
                    prior.append(struct.unpack('=d', byte)[0])
                    depth = []
                    res = []
                    for x in range(0,a):            
                        if x == 0:
                            depth.append(0)
                        byte = f.read(8)
                        b = struct.unpack('=d', byte)[0]
                        depth.append(b)
                    for x in range(0,a):
                        byte = f.read(8)
                        res.append(struct.unpack('=d', byte)[0])

                # store models ensembles
                px = np.ones([2*len(depth[1:]),2])
                px[0::2,0],px[1::2,0],px[1::2,1],px[2::2,1] = res[:],res[:],depth[1:],depth[1:-1]
                ensemble_models.append(px)


    finally:
        f.close()

    # save maximum likelihood
    maxLkh_mod = ensemble_models[lkh.index(max(lkh))]

    return ensemble_models, nLayers, lkh, maxLkh_mod, prior




def inputData(dat_file,appres=False):
    
    df_obs = pd.read_csv(dat_file, skiprows=0)
    
    C = 10000/(4*np.pi)

    nFreq = len(df_obs)
    data = np.zeros((nFreq, 4))
    data[:,0] = df_obs['freq'].values
    if appres:
        data[:,1] = np.exp(df_obs['Zr'].values)
        data[:,2] = np.exp(df_obs['Zi'].values)
    else:
        data[:,1] = C*np.exp(df_obs['Zr'].values)
        data[:,2] = C*np.exp(df_obs['Zi'].values)
    
    log_std = df_obs['std'].values
    # convert errors from log to linear scale (mean of errors in Zr and Zi)
    ZrERR = ((np.exp(df_obs['Zr'].values+log_std) - np.exp(df_obs['Zr'].values-log_std) )) / 2
    ZiERR = ((np.exp(df_obs['Zi'].values+log_std) - np.exp(df_obs['Zi'].values-log_std) )) / 2  
 
    data[:,3] = C * 0.5 * (ZrERR + ZiERR)

    return data




def sampResps(samp_resp_file, freq):
    
    C = 10000/(4*np.pi)
    nFreq = len(freq)
    f = open(samp_resp_file, "rb")
    stop = False
    ensemble_resps = []
    try:
        while True:
            Z_re = []
            Z_im = []
            for x in range(0,nFreq):
                byte = f.read(8)
                if not byte:
                    stop = True
                    break
                zr = struct.unpack('=d', byte)[0] * C
                Z_re.append(zr)
                
                byte = f.read(8)
                zi = struct.unpack('=d', byte)[0] * C
                Z_im.append(zi)
            if stop:
                break
            else:
                resps = np.zeros((nFreq,2))
                resps[:,0] = Z_re
                resps[:,1] = Z_im
                ensemble_resps.append(resps)
    finally:
        f.close()
    
    return ensemble_resps



"""
Conpute change point distributions 
Each function can be customized to filter out transitions

Input:
- single model from the ensemble 

Output:
- layers depth for this model

"""

@jit(nopython=False) 
def changePt(m):
    vals=[]
    #1 Discard first 10 meters
    #m1=m[m[:,1]>1]
    #1 Discard below ~5000 meters
    #m1=m1[m1[:,1]<3.7]
    # Loop over each layers:
    for lay in range(len(m)):
        vals.append(m[lay,1])
    return vals



@jit(nopython=False) 
def changePt_bsmt(m):
    vals=[]
    # Loop over each layers:
    for lay in range(len(m)-1):
        if m[lay,0] < 3 and m[lay+1,0] > 3:
            vals.append(m[lay,1])
    return vals




"""
Compute the 2D histograms to plot posterior model and response distributions 

Input:
- ensembles
- grid boundaries (resistivity/depth or data/frequency) 
- grid size

Output:
- 2D histogram
- stats (model only)
- chPts (model only)

"""

@jit(nopython=False) 
def histModels(ensemble_models,  grid_bounds, nx=300, nz=300):

    xmin = int(grid_bounds[0])
    xmax = int(grid_bounds[1])
    zmin = int(grid_bounds[2])
    zmax = int(grid_bounds[3])
    
    x_range = xmax-xmin
    hist = np.zeros((nz, nx))
    ZZ = np.linspace(zmin+(zmax-zmin)/(nz*2.), zmax+(zmax-zmin)/(nz*2.), nz+1)[:-1]
    models = np.zeros((nz,len(ensemble_models)))
    
    chPts=[]
    
    for i in range(0, len(ensemble_models)):
        p = 0
        m = ensemble_models[i][1::2] 
        vals0 = changePt(m) 
        chPts.append(vals0)

        for z in range(0,nz):
            while ZZ[z] > m[p,1]:
                p += 1
            v = int(math.ceil((((m[p,0]-xmin)/x_range)*nx)))
            hist[z,v-1] += 1
            models[z,i] = m[p,0]
    
    # calculate some statistics
    stats = np.zeros((nz,7))
    stats[:,0] = ZZ  

    for z in range(0,nz):
        hist[z,:] = hist[z,:] / sum(hist[z,:])

        stats[z,2] = np.median(models[z,:])
        stats[z,1] = np.percentile(models[z,:], 5)
        stats[z,3] = np.percentile(models[z,:], 95)
        stats[z,4] = np.std(models[z,:])
        stats[z,5] = np.mean(models[z,:])
        #stats[z,6] = float(scipy.stats.mode(models[z,:])[0])
    
    chPts = [item for sublist in chPts for item in sublist]

    return hist, stats, np.array(chPts)







    
@jit(nopython=False) 
def histResponses(ensemble_models, grid_bounds_resps, 
                  nRho=100, 
                  nT=100, 
                  nModels=1000,
                  datatype='Z'):    
    
    """
    Compute the posterior PDFfor the responses. To have a complete PDF even 
    in presence of gaps in the data we compute the fwd solution of nModels 
    to calculate the histogram (slow...)
    
    Input:
        - ensemble of models
        - grid boundaries:
                - [Tmin,Tmax,Rhomin,Rhomax,Phymin,Phymax] or
                     [Tmin,Tmax,Zrmin,Zrmax,Zimin,Zimax]
                - Y: frequency / period --> nT
                - X: resistivity in log10 (or impedance Z in log 10)) --> nRho
        - number of models to use for the plot
        - number of frequencies to compute the forward 
                (= as number of cells of the histograms, easier to code)
        
    Output:
        - normalized histograms for rho/phy or Zr/Zi  
                        (transposed for plotting)
    """
    
    if len(ensemble_models) < nModels:
        nModels = len(ensemble_models)
    randMod = random.sample(range(0, len(ensemble_models)), nModels)
    
    Tmin = grid_bounds_resps[0]
    Tmax = grid_bounds_resps[1]
    Rhomin = grid_bounds_resps[2]
    Rhomax = grid_bounds_resps[3]
    Phymin = grid_bounds_resps[4]
    Phymax = grid_bounds_resps[5]   
    Zmin = grid_bounds_resps[6]
    Zmax = grid_bounds_resps[7]  

    Rho_range = Rhomax-Rhomin
    Phy_range = Phymax-Phymin
    Z_range = Zmax-Zmin    

    histRho = np.zeros((nT, nRho))
    histPhy = np.zeros((nT, nRho))
    histZr = np.zeros((nT, nRho))
    histZi = np.zeros((nT, nRho))

    TT = np.linspace(Tmin+(Tmax-Tmin)/(nT*2.), Tmax+(Tmax-Tmin)/(nT*2.), nT+1)[:-1]
    C = 10000/(4*np.pi)
    
    for i in randMod:
        model = ensemble_models[i][1::2] 
        f, rho, phy, Z = MT.fwd1D(model,(1/10**TT))
        Zr = Z.real *C
        Zi = Z.imag *C
        
        for t in range(nT):
            if datatype == 'rhoPhy':
                v_rho = int(math.ceil((((np.log10(rho[t])-Rhomin)/Rho_range)*nRho)))
                histRho[t,v_rho-1] += 1
                v_phy = int(math.ceil((((phy[t]-Phymin)/Phy_range)*nRho)))
                histPhy[t,v_phy-1] += 1      
            elif datatype == 'Z':
                v_Zr = int(math.ceil((((np.log10(Zr[t])-Zmin)/Z_range)*nRho)))
                histZr[t,v_Zr-1] += 1    
                v_Zi = int(math.ceil((((np.log10(Zi[t])-Zmin)/Z_range)*nRho)))
                histZi[t,v_Zi-1] += 1    
    
    for T in range(0,nT):
        if datatype == 'rhoPhy':
            histRho[T,:] = histRho[T,:] / sum(histRho[T,:])    
            histPhy[T,:] = histPhy[T,:] / sum(histPhy[T,:])    
        elif datatype == 'Z':
            histZr[T,:] = histZr[T,:] / sum(histZr[T,:])     
            histZi[T,:] = histZi[T,:] / sum(histZi[T,:])     
    
    if datatype == 'rhoPhy':
        return np.flipud(histRho.T), np.flipud(histPhy.T)
    elif datatype == 'Z':
        return  np.flipud(histZr.T), np.flipud(histZi.T)




def read_invParams(params_file):
    with open(params_file) as f:
        lines = f.readlines()

    nChains = int(lines[1].split('=')[1][:-1])
    nIt = int(lines[2].split('=')[1][:-1])
    samples_perChain = int(lines[5].split('=')[1][:-1])
    rhoMin = int(lines[7].split('=')[1][:-1])
    rhoMax = int(lines[8].split('=')[1][:-1])
    
    return nChains, nIt, samples_perChain, rhoMin, rhoMax

# class inv(object):
#     def __init__(self,params_file):
#         with open(params_file) as f:
#             lines = f.readlines()
#         self.nChains = int(lines[1].split('=')[1][:-1])
#         self.nIt = int(lines[2].split('=')[1][:-1])
#         self.samples_perChain = int(lines[5].split('=')[1][:-1])
#         self.rhoMin = int(lines[7].split('=')[1][:-1])
#         self.rhoMax = int(lines[8].split('=')[1][:-1])
        

    