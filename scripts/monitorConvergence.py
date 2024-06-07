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
from os.path import exists
import sys
import numpy as np

sys.path.append("../src")
import ensembles


files_path = f'../projects/{project}/transdMT/outfolder'

site_ids = []
for file in os.listdir(f'{files_path}/csv'):
    if file.endswith(".csv") and not file.endswith("log.csv"): #and not file.startswith("AF3"):
        site_ids.append(file[:-4])
site_ids = np.sort(site_ids)


params_file = './inversionParameters.txt'
nChains, nIt, samples_perChain, rhoMin, rhoMax = ensembles.read_invParams(params_file)

print('Project: ',project)
for site_id in site_ids:             
    
    print(' Plotting convergeance statistics for MT site %s...'%site_id)
    
    thinning = 1000
    
    log_file = f'{files_path}/logs/{site_id}_log.csv'
    
    assert exists(log_file), '   convergence statistics for MT site %s not existing!!'%site_id
    df_obs = pd.read_csv(log_file, skiprows=0,header=None)
    for i in range(0,nChains):
        df_obs = df_obs.drop(i*nIt/thinning +i)
        
    df_obs = df_obs.reset_index()
    df_obs = np.array(df_obs, dtype=float)
    

    it_num = df_obs[:int(nIt/thinning),1]
    
    plt.figure(1,figsize=(8,4))
    
    plt.subplot(211)
    for i in range(0,nChains):
        # plt.semilogy(it_num,df_obs[' nll'][int(i*nIt/thinning):int((i+1)*nIt/thinning)],'-', lw=0.1)
        plt.semilogy(it_num,df_obs[int(i*nIt/thinning):int((i+1)*nIt/thinning),3],'-', lw=0.1)
    plt.xlim([0,nIt])
    plt.ylabel('Negative\nLog Likelihood')
    plt.title('MT site %s - MCMC Convergeance - %d chains - %dk it/chain'
              %(site_id, nChains, nIt/1000))
    plt.axvline(0.75*nIt, c='k', label = 'end of\nburn-in phase')
    plt.legend(loc=1, fontsize='small')

    plt.subplot(212)
    for i in range(0,nChains):
        # plt.semilogy(it_num,df_obs[' nll'][int(i*nIt/thinning):int((i+1)*nIt/thinning)],'-', lw=0.1)
        plt.plot(it_num,df_obs[int(i*nIt/thinning):int((i+1)*nIt/thinning),2],'-', lw=0.1)
    plt.xlim([0,nIt])
    plt.ylabel('Number of\ninterfaces')
    plt.xlabel('Iterations')
    plt.axvline(0.75*nIt, c='k')

    plt.tight_layout


    plt.savefig('%s/%s_convergence.png'%(files_path,site_id),dpi=300, bbox_inches="tight")
    plt.close('all')




