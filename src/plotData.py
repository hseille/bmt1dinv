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



import numpy as np
import matplotlib.pyplot as plt
import MT as MT


def plot_edi(site_id, dat1D, ss, dat, medfiltsp, 
        plot_rhoPhy=True, 
        plot_antidiag = True, 
        plot_diag = True,
        phase_90=False, 
        plot_size='medium',
        rho_lims = [-1, 4]):
    
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
    Z1DR = C*np.exp(dat1D['Z1DRLn'].values) 
    Z1DI = C*np.exp(dat1D['Z1DILn'].values) 
    ZrERR = ((np.exp(dat1D['Z1DRLn'].values+ss) - np.exp(dat1D['Z1DRLn'].values-ss) )) / 2
    ZiERR = ((np.exp(dat1D['Z1DILn'].values+ss) - np.exp(dat1D['Z1DILn'].values-ss) )) / 2  
    Zstd = C * 0.5 * (ZrERR + ZiERR)

    
    if plot_size == 'small':
        fig = plt.figure(4,figsize=(3, 6))
        ftsz =9
    elif plot_size == 'medium':
        fig = plt.figure(4,figsize=(5, 8))
        ftsz = 10
    
    if plot_rhoPhy:
        
        rhoXX, phyXX, drhoXX, dphyXX = MT.z2rhophy(f, Z['ZXXR'],Z['ZXXI'],Z['ZXX.VAR'])
        rhoXY, phyXY, drhoXY, dphyXY = MT.z2rhophy(f, Z['ZXYR'],Z['ZXYI'],Z['ZXY.VAR'])
        rhoYX, phyYX, drhoYX, dphyYX = MT.z2rhophy(f, Z['ZYXR'],Z['ZYXI'],Z['ZYX.VAR'])
        rhoYY, phyYY, drhoYY, dphyYY = MT.z2rhophy(f, Z['ZYYR'],Z['ZYYI'],Z['ZYY.VAR'])
        
        rho1D, phy1D, drho1D, dphy1D = MT.z2rhophy(f, Z1DR, Z1DI, (Zstd)**2)

        
        minimRho = np.log10(min(min(rhoXY),min(rhoYX)))
        maximRho = np.log10(max(max(rhoXY),max(rhoYX)))
        logERR_xx = ((np.log10(rhoXX+drhoXX)-np.log10(rhoXX-drhoXX)))*0.5
        logERR_xy = ((np.log10(rhoXY+drhoXY)-np.log10(rhoXY-drhoXY)))*0.5
        logERR_yx = ((np.log10(rhoYX+drhoYX)-np.log10(rhoYX-drhoYX)))*0.5
        logERR_yy = ((np.log10(rhoYY+drhoYY)-np.log10(rhoYY-drhoYY)))*0.5
        logERR_1D = ((np.log10(rho1D+drho1D)-np.log10(rho1D-drho1D)))*0.5
        
        plt.subplot(10,1,(1,3))
        if plot_antidiag:
            plt.errorbar(np.log10(1/f), np.log10(rhoXY), yerr = logERR_xy, 
                    fmt='r.',label= 'xy',zorder=32, 
                    elinewidth=0.6,markersize=5 ,
                    capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.2)
            plt.errorbar(np.log10(1/f), np.log10(rhoYX), yerr = logERR_yx, 
                    fmt='b.',label= 'yx',zorder=32, 
                    elinewidth=0.6,markersize=5 ,
                    capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.2)
        if plot_diag:
            plt.errorbar(np.log10(1/f), np.log10(rhoXX), yerr = logERR_xx, 
                    fmt='m.',label= 'xx',zorder=32, 
                    elinewidth=0.6,markersize=5 ,
                    capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.2)
            plt.errorbar(np.log10(1/f), np.log10(rhoYY), yerr = logERR_yy, 
                    fmt='g.',label= 'yy',zorder=32, 
                    elinewidth=0.6,markersize=5 ,
                    capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.2)
        plt.errorbar(np.log10(1/f), np.log10(rho1D), yerr = logERR_1D, 
                fmt='k.',label= '1D', zorder=31, 
                elinewidth=0.6,markersize=3 ,
                capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.8)
        plt.ylabel('Log$_{10}$ $\\rho_{app}$ ($\Omega$.m)',fontsize=ftsz)
        plt.title('%s'%(site_id),fontsize=ftsz)
        plt.grid(linestyle=':')
        plt.xticks(np.arange(-5,5,1),labels=[])
        plt.xlim(Tmin-0.5,Tmax+0.5)        
        plt.yticks(np.arange(-2,6,1),fontsize=ftsz)
        # plt.ylim(minimRho-0.5,maximRho+0.5)
        plt.ylim(rho_lims[0], rho_lims[1])
        plt.legend(fontsize=ftsz-2, loc=0, ncol=2)
        
        if phase_90:
            phyYX += 180
        
        plt.subplot(10,1,(4,6))
        if plot_antidiag:
            plt.errorbar(np.log10(1/f), phyXY, yerr = dphyXY, 
                    fmt='r.',zorder=32, 
                    elinewidth=0.6,markersize=5 ,
                    capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.2)
            plt.errorbar(np.log10(1/f), phyYX, yerr = dphyYX, 
                    fmt='b.',zorder=32, 
                    elinewidth=0.6,markersize=5 ,
                    capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.2)
        if plot_diag:
            plt.errorbar(np.log10(1/f), phyXX, yerr = dphyXX, 
                    fmt='m.',zorder=32, 
                    elinewidth=0.6,markersize=5 ,
                    capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.2)
            plt.errorbar(np.log10(1/f), phyYY, yerr = dphyYY, 
                    fmt='g.',zorder=32,
                    elinewidth=0.6,markersize=5 ,
                    capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.2)
        plt.errorbar(np.log10(1/f), phy1D, yerr = dphy1D, 
                fmt='k.',zorder=31, 
                elinewidth=0.6,markersize=3 ,
                capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.8)
        plt.ylabel('Phase ($^\circ$)',fontsize=ftsz)
        plt.grid(linestyle=':')
        plt.xticks(np.arange(-5,5,1),labels=[])
        plt.xlim(Tmin-0.5,Tmax+0.5)
        if phase_90:
            plt.yticks([0,45,90],fontsize=ftsz)
            plt.ylim(0,90);
        else:
            plt.yticks([-180,-90,0,90,180],fontsize=ftsz)
            plt.ylim(-180,180);
        
        plt.subplot(10,1,(7,8))
        #plt.plot(logT,beta, 'k.')
        plt.errorbar(np.log10(1/f), beta, yerr = 2*beta_err,
                fmt='k.', elinewidth=0.6,markersize=5 ,
                capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.8)
        plt.plot(np.log10(1/f), beta_filt, 'r-', lw=1,zorder=32, label="Median Filter")
        plt.ylabel('PT $\\beta$ ($^\circ$)',fontsize=ftsz)
        plt.grid(linestyle=':')
        plt.axhline(y=3, linestyle='-', color='grey')
        plt.axhline(y=-3, linestyle='-', color='grey')
        plt.xticks(np.arange(-5,5,1),labels=[])
        plt.xlim(Tmin-0.5,Tmax+0.5)
        plt.yticks([-10,-3,0,3,10],fontsize=ftsz)
        plt.ylim(-12, 12)
        plt.legend(loc=2,fontsize=ftsz-2, ncol=2)
        
        plt.subplot(10,1,(9,10))
        plt.errorbar(np.log10(1/f), ellip, yerr = 2*ellip_err, 
                fmt='k.', elinewidth=0.6,markersize=5 ,
                capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.8)
        plt.plot(np.log10(1/f), ellip_filt, 'r-', lw=1,zorder=32, label="Median Filter")
        plt.xlabel('Log$_{10}$ Period (sec)',fontsize=ftsz)
        plt.ylabel('PT ellip. $\\lambda$',fontsize=ftsz)
        plt.grid(linestyle=':')
        plt.axhline(y=0.1, linestyle='-', color='grey')
        plt.xticks(np.arange(-5,5,1),fontsize=ftsz)
        plt.xlim(Tmin-0.5,Tmax+0.5)
        plt.yticks([0,1],fontsize=ftsz)
        plt.ylim(-0.01, 1.2)
        plt.legend(loc=2,fontsize=ftsz-2, ncol=2)
    
    
    
    
    #IMPEDANCES 
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
        plt.errorbar(np.log10(1/f), beta, yerr = 2*beta_err,
                fmt='k.', elinewidth=0.6,markersize=5 ,
                capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.8)
        plt.plot(np.log10(1/f), beta_filt, 'r-', lw=1,zorder=32, label="Median Filter")
        plt.ylabel('PT $\\beta$ ($^\circ$)',fontsize=ftsz)
        plt.grid(linestyle=':')
        plt.axhline(y=3, linestyle='-', color='grey')
        plt.axhline(y=-3, linestyle='-', color='grey')
        plt.xticks(np.arange(-5,5,1),labels=[])
        plt.xlim(Tmin-0.5,Tmax+0.5)
        plt.yticks([-10,-3,0,3,10],fontsize=ftsz)
        plt.ylim(-12, 12)
        plt.legend(loc=2,fontsize=ftsz-2, ncol=2)
        
        plt.subplot(10,1,(9,10))
        plt.errorbar(np.log10(1/f), ellip, yerr = 2*ellip_err, 
                fmt='k.', elinewidth=0.6,markersize=5 ,
                capsize=2,capthick=0.6,mec='k',mew=0.5, alpha=0.8)
        plt.plot(np.log10(1/f), ellip_filt, 'r-', lw=1,zorder=32, label="Median Filter")
        plt.xlabel('Log$_{10}$ Period (sec)',fontsize=ftsz)
        plt.ylabel('PT ellip. $\\lambda$',fontsize=ftsz)
        plt.grid(linestyle=':')
        plt.axhline(y=0.1, linestyle='-', color='grey')
        plt.xticks(np.arange(-5,5,1),fontsize=ftsz)
        plt.xlim(Tmin-0.5,Tmax+0.5)
        plt.yticks([0,1],fontsize=ftsz)
        plt.ylim(-0.01, 1.2)
        plt.legend(loc=2,fontsize=ftsz-2, ncol=2)

    plt.tight_layout()
    
    
