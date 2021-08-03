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
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Hoël Seillé"
__email__ = "hoel.seille@csiro.au"
__status__ = "Beta"






# Read outputs of trans-d inversion


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
        
    #samp_models_file =r'C:\Users\sei029\Documents\Projects\BAY1DMT\Work/inversion_1D_transd\outputs\transdMT4\G05_tr9_EF10/133CsampModels'


    f = open(samp_models_file, "rb")
    
    lkh = []
    prior =[]
    nLayers = []
    # depths_of_change = []
    # positive_change = []
    # negative_change = []
    bb = 0
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
                    byte = f.read(8)                             # two lines added in version transdMT3 (prior H10)
                    prior.append(struct.unpack('=d', byte)[0])
                    depth = []
                    res = []
                    # d = 0
                    for x in range(0,a):            
                        if x == 0:
                            depth.append(0)
                        byte = f.read(8)
                        b = struct.unpack('=d', byte)[0]
                        #d = d + b                             # removed in version transdMT4 (we deal with depths now)
                        depth.append(b)
                    for x in range(0,a):
                        bbb = bb
                        byte = f.read(8)
                        bb = struct.unpack('=d', byte)[0]
                        res.append(struct.unpack('=d', byte)[0])
                #         if x != 0:
                #             if bb > bbb:
                #                 positive_change.append(1)
                #                 negative_change.append(0)
                #             else:
                #                 positive_change.append(0)
                #                 negative_change.append(-1)
                # for i in range(1, len(depth)-1):
                #     depths_of_change.append(depth[i])
    
                # SAVE IN ENSEMBLE OF MODELS
                
                px = np.ones([2*len(depth[1:]),2])
                px[0::2,0],px[1::2,0],px[1::2,1],px[2::2,1] = res[:],res[:],depth[1:],depth[1:-1]
                #plt.plot(px[:,0],-px[:,1],'k-',lw=1, label = 'True model')
                ensemble_models.append(px)

    
    finally:
        f.close()
        
    # layer_change = np.zeros((len(depths_of_change), 3))
    # layer_change[:,0] = depths_of_change
    # layer_change[:,1] = positive_change
    # layer_change[:,2] = negative_change
    
    # calculate maximum likelihood
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
        data[:,2] = np.exp(df_obs['Zr'].values)
    else:
        data[:,1] = C*np.exp(df_obs['Zr'].values)
        data[:,2] = C*np.exp(df_obs['Zi'].values)
    
    log_std = df_obs['std'].values
    # convert errors from log to linear scale (mean of errors in Zr and Zi... ?)
    ZrERR = ((np.exp(df_obs['Zr'].values+log_std) - np.exp(df_obs['Zr'].values-log_std) )) / 2
    ZiERR = ((np.exp(df_obs['Zi'].values+log_std) - np.exp(df_obs['Zi'].values-log_std) )) / 2  
 
    data[:,3] = C*0.5 * (ZrERR + ZiERR)

    return data

    
    
    

def sampResps(samp_resp_file, freq):
    
#work_dir = r'C:\Users\sei029\Documents\Projects\BAY1DMT\Work\inversion_1D_transd\outputs\test_plotting'
#siteid = 18
#samp_resp_file = '%s/transdin%ssampResps'%(work_dir,siteid)

    nFreq = len(freq)
    
    f = open(samp_resp_file, "rb")
    stop = 'False'
    ensemble_resps = []
    try:
        while True:
            Z_re = []
            Z_im = []
            for x in range(0,nFreq):
                byte = f.read(8)
                if not byte:
                    stop = 'True'
                    break
                zr = struct.unpack('=d', byte)[0]
                Z_re.append(zr)
                
                byte = f.read(8)
                zi = struct.unpack('=d', byte)[0]
                Z_im.append(zi)
            if stop == 'True':
                break
            else:
                a=1
                #ax1.semilogx(period, Z_re, color = 'cornflowerblue',lw=0.1,zorder=-32)
                #ax1.semilogx(period, Z_im, color = 'lightcoral',lw=0.1,zorder=-32)
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
Compute the posterior model and response distributions 

Input:
- ensembles
- grid boundaries (resistivity/depth or resistivity/frequency) 
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
    n = len(ensemble_models)
    hist = np.zeros((nz, nx))
    dz = zmax/nz
    ZZ = np.linspace(zmin+(zmax-zmin)/(nz*2.), zmax+(zmax-zmin)/(nz*2.), nz+1)[:-1]
    models = np.zeros((nz,len(ensemble_models)))
    
    chPts=[]
    #print('ok2')
    
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

# =============================================================================
#     plt.subplot(1,2,2)
#     plt.imshow(np.log10(hist), extent=[grid_bounds[0],grid_bounds[1],-grid_bounds[3],grid_bounds[2]], aspect='auto')
#     plt.plot(ensemble_models[i][:,0], -ensemble_models[i][:,1])
#     plt.ylim(-300,0)
#     plt.xlim(0.5,2)
# =============================================================================
    return hist, stats, np.array(chPts)#, changePt_hist







    
@jit(nopython=False) 
def histResponses(ensemble_models, grid_bounds_resps, nRho=100, nT=100, nModels=1000):    
    
    """
    Compute the posterior PDFfor the responses. To have a complete PDF even 
    in presence of gaps in the data we compute the fwd solution of nModels 
    to calculate the histogram (slow...)
    
    Input:
        - ensemble of models
        - grid boundaries ([Tmin,Tmax,Rhomin,Rhomax,Phymin,Phymax] / [Tmin,Tmax,Zrmin,Zrmax,Zimin,Zimax]])
                            (Y: frequency / period --> nT
                             X: resistivity in log10 (or impedance Z in log 10)) --> nRho
        - number of models to use for the plot
        - number of frequencies to compute the forward (= as number of cells of the histograms, easier to code)
        
    Output:
        - normalized histograms for rho/phy or Zr/Zi  (transposed for the plotting)
    """
    
    if len(ensemble_models) < nModels:
        nModels = len(ensemble_models)
    randMod = random.sample(range(0, len(ensemble_models)), nModels)
    #randMod = random.sample(range(0, len(ensemble_models)), 1)
    
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
    #n = nModels
    histRho = np.zeros((nT, nRho))
    histPhy = np.zeros((nT, nRho))
    histZr = np.zeros((nT, nRho))
    histZi = np.zeros((nT, nRho))
    #dz = zmax/nz
    TT = np.linspace(Tmin+(Tmax-Tmin)/(nT*2.), Tmax+(Tmax-Tmin)/(nT*2.), nT+1)[:-1]
    
    for i in randMod:
        model = ensemble_models[i][1::2] 
        f, rho, phy, Z = MT.fwd1D(model,(1/10**TT))
        
        for t in range(nT):
            v_rho = int(math.ceil((((np.log10(rho[t])-Rhomin)/Rho_range)*nRho)))
            histRho[t,v_rho-1] += 1
            v_phy = int(math.ceil((((phy[t]-Phymin)/Phy_range)*nRho)))
            histPhy[t,v_phy-1] += 1      
         
            v_Zr = int(math.ceil((((np.log10(Z[t]).real-Zmin)/Z_range)*nRho)))
            histZr[t,v_Zr-1] += 1    
            v_Zi = int(math.ceil((((np.log10(Z[t]).imag-Zmin)/Z_range)*nRho)))
            histZi[t,v_Zi-1] += 1    
    
    for T in range(0,nT):
        histRho[T,:] = histRho[T,:] / sum(histRho[T,:])    
        histPhy[T,:] = histPhy[T,:] / sum(histPhy[T,:])    
        histZr[T,:] = histZr[T,:] / sum(histZr[T,:])     
        histZi[T,:] = histZi[T,:] / sum(histZi[T,:])     
        
    return np.flipud(histRho.T),np.flipud(histPhy.T), np.flipud(histZr.T), np.flipud(histZi.T)


