'''
    File name: monitorConvergence.py
    Authors: Hoël Seillé / Gerhard Visser
    Date created: 25/01/2022
    Date last modified: 25/01/2022
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
# Define here the project
# =============================================================================

project = 'example'

# =============================================================================
# 
# =============================================================================



import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import numpy as np

sys.path.append("../src")
import ensembles


files_path = '../projects/%s/transdMT/outfolder'%project

site_ids = []
for file in os.listdir(files_path):
    if file.endswith(".csv") and not file.endswith("log.csv"):
        site_ids.append(file[:-4])


params_file = './inversionParameters.txt'
nChains, nIt, samples_perChain, rhoMin, rhoMax = ensembles.read_invParams(params_file)


for site_id in site_ids:             
    
    print('Plotting convergeance statistics for MT site %s...'%site_id)
    
    thinning = 1000
    
    log_file = '%s/%s_log.csv'%(files_path,site_id)
    df_obs = pd.read_csv(log_file, skiprows=0)
    for i in range(1,nChains):
        df_obs = df_obs.drop(i*nIt)
    df_obs = df_obs[::thinning]
    df_obs = df_obs.reset_index()
    df_obs = np.array(df_obs, dtype=float)
    
    #n_maxChains = 20
    it_num = np.linspace(0, nIt, int(nIt/thinning), dtype=int)
    
    plt.figure(1,figsize=(8,4))
    
    plt.subplot(211)
    for i in range(0,nChains):
        # plt.semilogy(it_num,df_obs[' nll'][int(i*nIt/thinning):int((i+1)*nIt/thinning)],'-', lw=0.1)
        plt.semilogy(it_num,df_obs[int(i*nIt/thinning):int((i+1)*nIt/thinning),2],'-', lw=0.1)
    plt.xlim([0,nIt])
    plt.ylabel('Normalized\nLog Likelihood')
    plt.title('MT site %s - MCMC Convergeance - %d chains - %dk it/chain'
              %(site_id, nChains, nIt/1000))
    plt.axvline(0.75*nIt, c='k', label = 'end of\nburn-in phase')
    plt.legend(loc=1, fontsize='small')

    plt.subplot(212)
    for i in range(0,nChains):
        # plt.semilogy(it_num,df_obs[' nll'][int(i*nIt/thinning):int((i+1)*nIt/thinning)],'-', lw=0.1)
        plt.plot(it_num,df_obs[int(i*nIt/thinning):int((i+1)*nIt/thinning),1],'-', lw=0.1)
    plt.xlim([0,nIt])
    plt.ylabel('Number of\ninterfaces')
    plt.xlabel('Iterations')
    plt.axvline(0.75*nIt, c='k')

    plt.tight_layout


    plt.savefig('%s/%s_convergeance.png'%(files_path,site_id),dpi=300, bbox_inches="tight")
    plt.close('all')











# fig = plt.figure(10,figsize=(10,7),constrained_layout=True)
# from matplotlib.gridspec import GridSpec
# gs = GridSpec(3, 6, figure=fig)
# ax1 = fig.add_subplot(gs[0, :4])
# ax2 = fig.add_subplot(gs[0, 4:])
# ax3 = fig.add_subplot(gs[1, :4])
# ax4 = fig.add_subplot(gs[1, 4:])
# ax5 = fig.add_subplot(gs[2, :4])
# ax6 = fig.add_subplot(gs[2, 4:])


# rms=[]
# for mod in range(len(ensemble_resps)):
#     Zr = ensemble_resps[mod][:,0]
#     Zi = ensemble_resps[mod][:,1]
#     resp = np.array([Zr.T,Zi.T]).T
#     rms.append(MT.rms(obs_dat,resp))


# for i in range(1,n_maxChains):
#     for j in range(1,nChains):
#         ax1.plot(it_num, rms[i*samples_perChain:(i+1)*samples_perChain], '-', lw=0.1)
# ax1.set_xlim([0.75 * nIt, nIt])
# # plt.ylim([df_obs[' nll'][750000:].mean()-5*df_obs[' nll'][750000:].std(),df_obs[' nll'][750000:].mean()+10*df_obs[' nll'][750000:].std()])
# ax1.set_ylabel('RMS\n ')
# ax1.set_title('%s'%site_id)


# weights = np.ones_like(rms) / len(rms)
# ax2.hist(rms,20,
#          color='k', 
#          histtype='bar', 
#          ec='white',
#          align='left', 
#          weights=weights)
# ax2.set_xlabel('RMS')
# ax2.set_ylabel('probability')

# for i in range(1,n_maxChains):
#     for j in range(1,nChains):
#         ax3.plot(it_num, lkh[i*samples_perChain:(i+1)*samples_perChain], '-', lw=0.1)
# ax3.set_xlim([0.75 * nIt, nIt])
# # plt.ylim([df_obs[' nll'][750000:].mean()-5*df_obs[' nll'][750000:].std(),df_obs[' nll'][750000:].mean()+10*df_obs[' nll'][750000:].std()])
# ax3.set_ylabel('Normalized\nLog Likelihood')

# weights = np.ones_like(lkh) / len(lkh)
# ax4.hist(lkh,len(np.arange(np.array(lkh).min()-1, np.array(lkh).max()+1, 1))-2,
#          color='k', 
#          histtype='bar', 
#          ec='white',
#          align='left', 
#          weights=weights)
# ax4.set_xlabel('Normalized Log Likelihood')
# ax4.set_ylabel('probability')

# for i in range(1,n_maxChains):
#     for j in range(1,nChains):
#         ax5.plot(it_num, nLayers[i*samples_perChain:(i+1)*samples_perChain], '-', lw=0.1)
# ax5.set_xlim([0.75 * nIt, nIt])
# # plt.ylim([0,df_obs['layers'][750000:].mean()+5*df_obs['layers'][750000:].std()])
# ax5.set_ylabel('Number of\ninterfaces')
# ax5.set_xlabel('Iterations')

# weights = np.ones_like(nLayers) / len(nLayers)
# ax6.hist(nLayers,len(np.arange(np.array(nLayers).min()-1, np.array(nLayers).max()+1, 1))-2,
#          color='k', 
#          histtype='bar', 
#          ec='white',
#          align='left', 
#          weights=weights)
# ax6.set_xlabel('# interfaces')
# ax6.set_ylabel('probability')


# plt.savefig('%s/%s_stats.png'%(files_path,site_id),dpi=300, bbox_inches="tight")
# plt.close('all')
