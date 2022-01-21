'''
    File name: plotPDFs.py
    Authors: Hoël Seillé / Gerhard Visser
    Date created: 01/10/2020
    Date last modified: 20/01/2022
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
import matplotlib.pyplot as plt
import matplotlib




def models(siteid,hist_norm_model,grid_bounds,stats, chPts):
    
    
    """
    Plot PDF of models

    Input:
    All the inputs are outputs of the functions posteriorHistograms.py
    - site Id
    - histogram
    - grid_bounds: boundaries of the histogram (min/max resistivity and depth)
    - ensemble stats 
    
    Output:
    - plot
    
    Options:
    - save plot ? no
    
    """
    fig = plt.figure(1, figsize=(4,4.2))


    ax1 = fig.add_axes([0.18,0.1,0.50,0.8])   # PDF model
    ax2 = fig.add_axes([0.7,0.1,0.15,0.8])   # chpT
    ax3 = fig.add_axes([0.87,0.1,0.02,0.8])   # colorbar
    
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])

    # ax1.text(-4.5,30, 'b)', fontsize=12)

    color = 'jet'
    jetBig = matplotlib.cm.get_cmap(color, 256)
    new_color = jetBig(np.linspace(0.05, 0.85, 256))
    endTransparency=np.linspace(0,1,100)
    new_color[:100,3] = endTransparency
    color = matplotlib.colors.ListedColormap(new_color)
  
    # Clip histogram for better visualisation
    h1 = np.clip(hist_norm_model, 0.001, 0.2)

    im = ax1.imshow(np.log10(h1), 
                    interpolation = None, 
                    cmap=color, 
                    extent=[int(grid_bounds[0]),int(grid_bounds[1]),
                            -int(grid_bounds[3]),-int(grid_bounds[2])], 
                    aspect='auto', 
                    alpha=0.75)
    
    ax1.plot(stats[:,2], -stats[:,0],'k',
             lw=1, label = 'Median', zorder = 32)
    
    ax1.plot(stats[:,1], -stats[:,0],'k--',
             lw=1, label = '5th/95th perc.', zorder = 31)
    
    ax1.plot(stats[:,3], -stats[:,0],'k--',
             lw=1, zorder = 31)

    # ax1.set_xlim(1,5)
    ax1.set_ylim(-int(grid_bounds[3]),0)
    ax1.grid(linestyle=':')
    ax1.set_ylabel('Depth (m)');
    ax1.set_xlabel('Log$_{10}$ Resistivity ($\Omega$.m)')
    # ax1.tick_params(axis='x',top='on')
    # ax1.tick_params(axis='x',labeltop='off')
    ax1.set_title('%s'%siteid)
    #ax1.legend(loc=3,framealpha = 1)
    
    nbins = int(abs(int(grid_bounds[3])-int(grid_bounds[2])) / 100)
    
    bins = np.arange(-int(grid_bounds[3]),-int(grid_bounds[2])+1,nbins)
    ax2.grid(linestyle=':')
    ax2.hist(-chPts, bins=bins,  orientation="horizontal", 
                                 density=False,
                                 alpha=1, histtype='bar', 
                                 ec='k',color='whitesmoke')
    ax2.set_ylim(-int(grid_bounds[3]),0)
    ax2.set_xlabel('\nDensity')
    ax2.set_title('Change\nPoints',fontsize=10)
    
    fig.colorbar(im,cax=ax3, 
                 ticks=[0,-1,-2,-3,-4],
                 label = 'Log$_{10}$ PDF')#, orientation = 'horizontal')

    #plt.tight_layout()
    
   
    
    
    

def responses(siteid, histResp, grid_bounds_resps, datatype='rho',colorbar=True):
    
    """
    Plot PDF of responses

    Input:
    All the inputs are outputs of the functions posteriorHistograms.py
    - site Id
    - histogram
    - grid_bounds: boundaries of the histogram (min/max resistivity and depth)

    
    Output:
    - plot
    
    Options:
    - save plot ? no
    
    """
    #histResp = histPhy
    #datatype='phy'
    #fig = plt.figure(10)
    
    color = 'jet'
    jetBig = matplotlib.cm.get_cmap(color, 256)
    new_color = jetBig(np.linspace(0.05, 0.85, 256))
    endTransparency=np.linspace(0,1,100)
    new_color[:100,3] = endTransparency
    color = matplotlib.colors.ListedColormap(new_color)
  
    # Clip histogram for better visualisation
    h1 = np.clip(histResp, 0.001, 0.2)
    
    if datatype == 'rho':  
        grid_bounds = [grid_bounds_resps[0],grid_bounds_resps[1],
                       grid_bounds_resps[2],grid_bounds_resps[3]]
    if datatype == 'phy':  
        grid_bounds = [grid_bounds_resps[0],grid_bounds_resps[1],
                       grid_bounds_resps[4],grid_bounds_resps[5]]
    if datatype == 'Z':  
        grid_bounds = [grid_bounds_resps[0],grid_bounds_resps[1],
                       grid_bounds_resps[6],grid_bounds_resps[7]]

    im = plt.imshow(np.log10(h1), interpolation = 'bicubic', 
                    cmap=color, 
                    extent=[grid_bounds[0],grid_bounds[1],
                            grid_bounds[2],grid_bounds[3]], 
                    aspect='auto')
    
    if datatype == 'rho':  
        plt.ylabel('Log$_{10}$ App. Res.($\Omega$.m)');plt.title('%s'%siteid)
    if datatype == 'phy':  
        plt.ylabel('Phase (degrees)'), plt.yticks([0,45,90])
    if datatype == 'Z':    
        plt.ylabel('Log$_{10}$ Z ($\Omega$)');plt.title('MT%s'%siteid)
    #plt.xlabel('Log$_{10}$ Period (sec)')
    plt.xticks(np.arange(-5,5))
    
    plt.grid(linestyle=':')
    #plt.xlim=(grid_bounds_resps[0],grid_bounds_resps[1])
   

    #plt.legend(loc=3,framealpha = 1)
    plt.colorbar(im,ticks=[0,-1,-2,-3,-4],label = 'Log$_{10}$ PDF')
    #fig.colorbar(im,ticks=[0,-1,-2,-3,-4])
    
    #plt.tight_layout()    
    
    
