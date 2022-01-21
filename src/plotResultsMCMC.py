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





# =============================================================================
# Define here the parameters and folders to create the results plots
# =============================================================================

project = 'example'
DepthMin = 0
DepthMax = 1500
PlotModels = True
PlotResponses = True
plot_niblettBostick =False
DepthLogScale = False  # plt the depth in log scale
files_path = '../projects/%s/transdMT/outfolder'%project

# =============================================================================
# 
# =============================================================================










import sys
import numpy as np
import matplotlib.pyplot as plt
import os
import copy

import ensembles as ensembles
import MT as MT
import plotPDFs as plotPDFs


print('Plot models: ',PlotModels)
print('Plot responses: ',PlotResponses)
print('Plot Niblett-Bostick depth-transform: ',plot_niblettBostick)


site_ids = []
for file in os.listdir(files_path):
    if file.endswith(".csv"):
        site_ids.append(file[:-4])


for site_id in site_ids:
    
    if PlotModels:
        
        samp_models_file =r'%s/%sCsampModels'%(files_path,site_id)
        
        print('\nMT site %s...'%site_id)
        
        print('    loading models ensemble...')
        # load ensemble of models
        ensemble_models, nLayers, lkh, maxLkh_mod, prior  = ensembles.sampModels(samp_models_file)
        #define boundaries of the histogram
        if DepthLogScale:  DepthMax = np.log10(DepthMax); DepthMin=1
        grid_bounds = [-2,6,DepthMin,DepthMax]
        
        #transform resistivities to log10
        LogEnsemble_models = copy.deepcopy(ensemble_models)
        for i in range(len(ensemble_models)):
            LogEnsemble_models[i][:,0] = np.log10(ensemble_models[i][:,0])
            if DepthLogScale:  LogEnsemble_models[i][:,1] = np.log10(ensemble_models[i][:,1])
            
        print('    computing and plotting model PDF...')
        #compute histograms of the posteriors for the models
        hist_norm_model, stats, chPts = ensembles.histModels(LogEnsemble_models,  
                                                             grid_bounds, 
                                                             nx=100, 
                                                             nz=300)
        #plot PDF model
        plotPDFs.models(site_id,hist_norm_model,grid_bounds,stats, chPts[chPts>10])
        # if DepthLogScale: plt.ylabel('Log$_{10}$ Depth (m)')
        # else: plt.ylabel('Depth (m)')
        plt.savefig('%s/%s_model.png'%(files_path,site_id),dpi=300, bbox_inches="tight")
        plt.close('all')


    
    
    
    if PlotResponses:
        #==============================================================================
        # # PDF RESPONSES
        #==============================================================================
        
        # read input data
        dat_file = r'%s/%s.csv'%(files_path,site_id)
        obs_dat = ensembles.inputData(dat_file, appres=False)
        ff = obs_dat[:,0]
        logT= np.log10(1/ff)
        nfreq=len(logT)
        
        
        rho, phy, drho, dphy = MT.z2rhophy(obs_dat[:,0],obs_dat[:,1],
                                           obs_dat[:,2],obs_dat[:,3]**2)

        samp_resps_file = r'%s/%sCsampResps'%(files_path,site_id)
        ensemble_resps = ensembles.sampResps(samp_resps_file, ff)
        
        #grid boundaries ([Tmin,Tmax,Rhomin,Rhomax,Phymin,Phymax,Zmin,Zmax] 
        Tmin=logT.min()
        Tmax=logT.max()
        grid_bounds_resps = [Tmin-0.0,Tmax+0.0,0,4,0,90,-5,2]    
        
        print('    computing and plotting response PDF...')
        #compute histograms of the posteriors for the responses
        histRho,histPhy, histZr,histZi = ensembles.histResponses(ensemble_models, 
                                                                 grid_bounds_resps,
                                                                 nRho=300, 
                                                                 nT=100, 
                                                                 nModels=100)
        
        # NORMAL DATA
        fig = plt.figure(2, figsize=(4,5))
        #plot PDF model
        plt.subplot(7,1,(1,3))
        plotPDFs.responses(site_id,histRho,grid_bounds_resps, datatype = 'rho')
        logERR = ((np.log10(rho+drho)-np.log10(rho-drho)))*0.5
        plt.errorbar(logT, np.log10(rho), yerr = logERR, 
                     fmt='k+',zorder=32, elinewidth=0.7,
                     markersize=2, capsize=0)

        plt.subplot(7,1,(4,5))
        plotPDFs.responses(site_id,histPhy,grid_bounds_resps, datatype = 'phy')
        plt.errorbar(logT, phy, yerr = dphy, 
                     fmt='k+',zorder=32, elinewidth=0.7,
                     markersize=2, capsize=0)

        
        if plot_niblettBostick:
            ax=plt.subplot(7,1,(6,7))
            rho_nb, depth_nb = MT.niblettBostick_depthTransform(rho, phy, 1/ff)
            #plt.loglog(depth_nb, rho_nb,'ko')
            plt.semilogy(logT, depth_nb,'k.')
            #plt.plot( np.log10(rho_nb),-depth_nb,'ko')
            plt.grid(linestyle=':')
            plt.ylabel('Depth (m)');plt.ylim([1, 2*depth_nb[-1]])
            plt.xlabel('Log$_{10}$ Period (sec)')
            plt.yticks([10,100,1000,10000,100000])
            plt.colorbar(ticks=[],drawedges=False, shrink=0.001)
            plt.title('Niblett-Bostick Depth Transform',fontsize=12)
            
            plt.ylim(10,100000)
            plt.xticks(np.arange(-5,5))
            plt.xlim(grid_bounds_resps[0],grid_bounds_resps[1])
        else: 
            plt.xlabel('Log$_{10}$ Period (sec)')
            
        
        plt.tight_layout() 
        plt.savefig('%s/%s_fit.png'%(files_path,site_id),dpi=300, bbox_inches="tight")
        plt.close('all')

