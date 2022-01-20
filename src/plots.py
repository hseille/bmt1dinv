'''
    File name: plots.py
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



#Plotting tools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MT as MT


def plot_edi(site_id, dat1D, ss, dat, medfiltsp, plot_rho=True):
    
    Z = dat[0]
    beta = dat[1][:,4]
    ellip = dat[1][:,5]
    beta_err = dat[2][:,3]
    ellip_err = dat[2][:,4]
    
    f = Z['FREQ'].values
    Tmin = min(np.log10(1 / f))
    Tmax = max(np.log10(1 / f))
    
    ellip_filt = MT.medianFilter(f, ellip, medfiltsp)
    beta_filt = MT.medianFilter(f, beta, medfiltsp)
    
    C = 10000/(4*np.pi)
    ZdetR = C*np.exp(dat1D['ZdetRLn'].values) 
    ZdetI = C*np.exp(dat1D['ZdetILn'].values) 
    ZrERR = ((np.exp(dat1D['ZdetRLn'].values+ss) - np.exp(dat1D['ZdetRLn'].values-ss) )) / 2
    ZiERR = ((np.exp(dat1D['ZdetILn'].values+ss) - np.exp(dat1D['ZdetILn'].values-ss) )) / 2  
    Zstd = C * 0.5 * (ZrERR + ZiERR)
    # print(Zstd)
    
    fig = plt.figure(4,figsize=(5, 8))
    
    if plot_rho:
        
        rhoXX, phyXX, drhoXX, dphyXX = MT.z2rhophy(f, Z['ZXXR'],Z['ZXXI'],Z['ZXX.VAR'])
        rhoXY, phyXY, drhoXY, dphyXY = MT.z2rhophy(f, Z['ZXYR'],Z['ZXYI'],Z['ZXY.VAR'])
        rhoYX, phyYX, drhoYX, dphyYX = MT.z2rhophy(f, Z['ZYXR'],Z['ZYXI'],Z['ZYX.VAR'])
        rhoYY, phyYY, drhoYY, dphyYY = MT.z2rhophy(f, Z['ZYYR'],Z['ZYYI'],Z['ZYY.VAR'])
        
        rho1D, phy1D, drho1D, dphy1D = MT.z2rhophy(f, ZdetR, ZdetI, (Zstd)**2)
        
        # print(drho1D**0.5)
        

        
        minimRho = np.log10(min(min(rhoXY),min(rhoYX)))
        maximRho = np.log10(max(max(rhoXY),max(rhoYX)))
        logERR_xx = ((np.log10(rhoXX+drhoXX)-np.log10(rhoXX-drhoXX)))*0.5
        logERR_xy = ((np.log10(rhoXY+drhoXY)-np.log10(rhoXY-drhoXY)))*0.5
        logERR_yx = ((np.log10(rhoYX+drhoYX)-np.log10(rhoYX-drhoYX)))*0.5
        logERR_yy = ((np.log10(rhoYY+drhoYY)-np.log10(rhoYY-drhoYY)))*0.5
        logERR_1D = ((np.log10(rho1D+drho1D)-np.log10(rho1D-drho1D)))*0.5
        
        plt.subplot(10,1,(1,3))
        plt.errorbar(np.log10(1/f), np.log10(rhoXX), 
                     yerr = logERR_xx, fmt='m.',label= 'xx',zorder=32, 
                     elinewidth=1,markersize=5, alpha=0.3)
        plt.errorbar(np.log10(1/f), np.log10(rhoXY), 
                     yerr = logERR_xy, fmt='r.',label= 'xy',zorder=32, 
                     elinewidth=1,markersize=5)
        plt.errorbar(np.log10(1/f), np.log10(rhoYX), 
                     yerr = logERR_yx, fmt='b.',label= 'yx',zorder=32, 
                     elinewidth=1,markersize=5)
        plt.errorbar(np.log10(1/f), np.log10(rhoYY), 
                     yerr = logERR_yy, fmt='g.',label= 'yy',zorder=32, 
                     elinewidth=1,markersize=5, alpha=0.3)
        plt.errorbar(np.log10(1/f), np.log10(rho1D), 
                     yerr = logERR_1D, fmt='k.',label= 'det',zorder=1, 
                     elinewidth=1,markersize=5, capsize=0,alpha=0.5)
        plt.ylabel('Log$_{10}$ $\\rho_{app}$ ($\Omega$.m)')
        plt.title('%s'%(site_id),fontsize=10)
        #plt.tick_params(axis='x',labelbottom='off')
        plt.grid(linestyle=':')
        plt.xlim(Tmin-0.5,Tmax+0.5)
        plt.yticks(np.arange(-2,6,1))
        plt.ylim(minimRho-0.5,maximRho+0.5)
        plt.legend(fontsize=8, loc=0)
        
        
        plt.subplot(10,1,(4,6))
        plt.errorbar(np.log10(1/f), phyXX, 
                     yerr = dphyXX, fmt='m.',zorder=32, 
                     elinewidth=1,markersize=5,  alpha=0.3)
        plt.errorbar(np.log10(1/f), phyXY, 
                     yerr = dphyXY, fmt='r.',zorder=32, 
                     elinewidth=1,markersize=5)
        plt.errorbar(np.log10(1/f), phyYX, 
                     yerr = dphyYX, fmt='b.',zorder=32, 
                     elinewidth=1,markersize=5, )
        plt.errorbar(np.log10(1/f), phyYY, 
                     yerr = dphyYY, fmt='g.',zorder=32, 
                     elinewidth=1,markersize=5, alpha=0.3)
        plt.errorbar(np.log10(1/f), phy1D, 
                     yerr = dphy1D, fmt='k.',zorder=1, 
                     elinewidth=1,markersize=5, capsize=0,alpha=0.5)
        plt.ylabel('Phase (degrees)')
        #plt.title('PHASE',fontsize=10)
        plt.grid(linestyle=':')
        plt.xlim(Tmin-0.5,Tmax+0.5)
        plt.yticks([-180,-90,0,90,180])
        plt.ylim(-180,180);
        #plt.legend(fontsize=8)
        
        plt.subplot(10,1,(7,8))
        #plt.plot(logT,beta, 'k.')
        plt.errorbar(np.log10(1/f), beta, 
                     yerr = 2*beta_err, fmt='k.', elinewidth=1,markersize=7)
        plt.plot(np.log10(1/f), beta_filt, 
                 'r-', zorder=32, label="Median Filter")
        plt.ylabel('$\\beta$ (degrees)')
        plt.grid(linestyle=':')
        plt.axhline(y=3, linestyle='-', color='grey')
        plt.axhline(y=-3, linestyle='-', color='grey')
        plt.xlim(Tmin-0.5,Tmax+0.5)
        plt.ylim(-12, 12)
        plt.legend(loc=2,fontsize=8)
        plt.title('PHASE TENSOR $\\beta$',fontsize=10)
        
        plt.subplot(10,1,(9,10))
        plt.errorbar(np.log10(1/f), ellip, 
                     yerr = 2*ellip_err, fmt='k.', elinewidth=1,markersize=7)
        plt.plot(np.log10(1/f), ellip_filt, 
                 'r-', zorder=32, label="Median Filter")
        plt.xlabel('Log$_{10}$ Period (sec)')
        plt.ylabel('$\\lambda$')
        plt.grid(linestyle=':')
        plt.axhline(y=0.1, linestyle='-', color='grey')
        plt.xlim(Tmin-0.5,Tmax+0.5)
        plt.ylim(-0.01, 1.2)
        plt.legend(loc=2,fontsize=8)
        plt.title('ELLIPTICITY  $\\lambda$',fontsize=10)
    
    
    
    
    #TO PLOT IMPEDANCES!!!  uncomment
    else:
        plt.subplot(10,1,(1,3))
        plt.errorbar(np.log10(1/f), np.log(Z['ZXYR'].values), 
                     yerr = np.log(1+((1*((Z['ZXY.VAR'])**0.5))/Z['ZXYR'].values)), 
                     fmt='r.',label= 'Zxy',zorder=32, elinewidth=1,markersize=7)
        plt.errorbar(np.log10(1/f), np.log(-Z['ZYXR'].values), 
                     yerr =  np.log(1+((1*((-Z['ZXY.VAR'])**0.5))/Z['ZYXR'].values)), 
                     fmt='b.',label= 'Zyx',zorder=32, elinewidth=1,markersize=7)
        plt.ylabel('Log Zr ($\Omega$)')
        plt.title('IMPEDANCE Z (Re)',fontsize=10)
        #plt.tick_params(axis='x',labelbottom='off')
        plt.grid(linestyle=':')
        plt.xlim(Tmin-0.5,Tmax+0.5)
        plt.legend(fontsize=8)
        
        plt.subplot(10,1,(4,6))
        plt.errorbar(np.log10(1/f), np.log(Z['ZXYI'].values), 
                     yerr =  np.log(1+((1*((Z['ZXY.VAR'])**0.5))/Z['ZXYI'].values)), 
                     fmt='r.',label= 'Zxy',zorder=32, elinewidth=1,markersize=7)
        plt.errorbar(np.log10(1/f), np.log(-Z['ZYXI'].values), 
                     yerr =  np.log(1+((1*((-Z['ZXY.VAR'])**0.5))/Z['ZYXI'].values)), 
                     fmt='b.',label= 'Zyx',zorder=32, elinewidth=1,markersize=7)
        
        plt.ylabel('Log -Zi ($\Omega$)')
        plt.title('IMPEDANCE Z (Im)',fontsize=10)
        plt.grid(linestyle=':')
        plt.xlim(Tmin-0.5,Tmax+0.5)
        plt.legend(fontsize=8)
        
        plt.subplot(10,1,(7,8))
        #plt.plot(logT,beta, 'k.')
        plt.errorbar(np.log10(1/f), beta, 
                     yerr = 2*beta_err, fmt='k.', elinewidth=1,markersize=7)
        plt.plot(np.log10(1/f), beta_filt, 
                 'r-', zorder=32, label="Median Filter")
        plt.ylabel('$\\beta$ (degrees)')
        plt.grid(linestyle=':')
        plt.axhline(y=3, linestyle='-', color='grey')
        plt.axhline(y=-3, linestyle='-', color='grey')
        plt.xlim(Tmin-0.5,Tmax+0.5)
        plt.ylim(-12, 12)
        plt.legend(loc=2,fontsize=8)
        plt.title('PHASE TENSOR $\\beta$',fontsize=10)
        
        plt.subplot(10,1,(9,10))
        plt.errorbar(np.log10(1/f), ellip, 
                     yerr = 2*ellip_err, fmt='k.', elinewidth=1,markersize=7)
        plt.plot(np.log10(1/f), ellip_filt, 
                 'r-', zorder=32, label="Median Filter")
        plt.xlabel('Log$_{10}$ Period (sec)')
        plt.ylabel('$\\lambda$')
        plt.grid(linestyle=':')
        plt.axhline(y=0.1, linestyle='-', color='grey')
        plt.xlim(Tmin-0.5,Tmax+0.5)
        plt.ylim(-0.01, 1.2)
        plt.legend(loc=2,fontsize=8)
        plt.title('ELLIPTICITY  $\\lambda$',fontsize=10)

    plt.tight_layout()
    
    
