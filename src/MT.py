'''
    File name: MT.py
    Authors: Hoël Seillé / Gerhard Visser
    Date created: 01/06/2020
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


import numpy as np
import pandas as pd




def readEDI(file_name):
    
    """
    function that reads EDI files
    Z in linear scale; format: ohm
    
    Input:
    - complete path to the EDI file
    
    Output:
    - data (EDI pandas DF format)
    - site_id (string)
    - coord (panda DF of 7 values (['Lat_deg']['Lat_min']['Lat_sec']['Long_deg']['Long_min']['Long_sec']['Elev'])
             
    Pandas DF format:        
    # ['FREQ',
    # 'ZXXR','ZXXI','ZXX.VAR',
    # 'ZXYR','ZXYI','ZXY.VAR',
    # 'ZYXR','ZYXI','ZYX.VAR',
    # 'ZYYR','ZYYI','ZYY.VAR',
    # 'TXR.EXP','TXI.EXP','TXVAR.EXP',
    # 'TYR.EXP','TYI.EXP','TYVAR.EXP']
    
    """
    
    import numpy as np
    import re
    import pandas as pd

    coord = pd.DataFrame({'Lat_deg':0,'Lat_min':0,'Lat_sec':0,'Long_deg':0,'Long_min':0,'Long_sec':0,'Elev':0},index=[0])

    with open(file_name, 'r') as f:
        data = f.readlines()
        for i,line in enumerate(data):
            line = data[i]
            words = line.split()
            
            #READ SITE NAME
            if any("DATAID" in s for s in words):
                words = ''.join(words)
                #print words
                site_id = re.search('\"([^"]+)', words).group(1)

            #READ NUMBER OF FREQUENCIES
            if any("NFREQ" in s for s in words):
                words = ''.join(words)
                nfreq_str = (re.findall('\d+', words))
                nfreq = int(nfreq_str[0])

             #READ LATITUDE
            if any("REFLAT" in s for s in words):
                words = ''.join(words)
                reflat_str = (re.findall('\-?\d+', words))
                coord['Lat_deg'] = str(reflat_str[0])
                coord['Lat_min'] = str(reflat_str[1])
                coord['Lat_sec'] = str('.'.join(reflat_str[2:]))				

            #READ LONGITUDE
            if any("REFLONG" in s for s in words):
                words = ''.join(words)
                reflong_str = (re.findall('\-?\d+', words))
                coord['Long_deg'] = str(reflong_str[0])
                coord['Long_min'] = str(reflong_str[1])
                coord['Long_sec'] = str('.'.join(reflong_str[2:]))	

            #READ ELEVATION
            if any("REFELEV" in s for s in words):
                words = ''.join(words)
                refelev_str = (re.findall('\-?\d+', words))
                coord['Elev'] = str('.'.join(refelev_str[0:]))				

    #READ MT DATA
    param = ['FREQ',
             'ZXXR','ZXXI','ZXX.VAR',
             'ZXYR','ZXYI','ZXY.VAR',
             'ZYXR','ZYXI','ZYX.VAR',
             'ZYYR','ZYYI','ZYY.VAR',
             'TXR.EXP','TXI.EXP','TXVAR.EXP',
             'TYR.EXP','TYI.EXP','TYVAR.EXP',]
    
    edi_data = np.empty((nfreq, len(param)))
    
    with open(file_name, 'r') as f:
        data = f.readlines()
        for i,line in enumerate(data):
            line = data[i]
            words = line.split()
            
            for col, data_type in enumerate(param):
                aa=[]            
                if ('>%s' %data_type) in words:
                    for k in range (1,1000):                   
                        if any(">" in s for s in data[i+k].split()):
                                break
                        else:
                            a = data[i+k].split()
                            aa += a
                    edi_data[:,col] = aa
    
    
    # write to Pandas format
    edi_pd = pd.DataFrame(edi_data)
    edi_pd.columns = param

    return (edi_pd, site_id, coord)





def phaseTensor(Z):
    
    # MAKE THAT FUNCTION MORE GENERAL AND SHORTER
    """
    Function that calculates phase tensor and parameters
    Following Caldwell et al. (2004 GJI) / Bibby et al. (2005 GJI)
    
    Input:
    - data (pandas DF format for Z)
    
    Output:
    - phase tensor matrix
    - phase tensor parameters (phmax,phmin,alpha,beta,ellip,azimuth)
    
    """   

    if isinstance(Z, pd.DataFrame):
        freq = Z['FREQ'].values
        nf = len(freq)
        X11 = Z['ZXXR'].values
        Y11 = Z['ZXXI'].values
        X12 = Z['ZXYR'].values
        Y12 = Z['ZXYI'].values
        X21 = Z['ZYXR'].values
        Y21 = Z['ZYXI'].values
        X22 = Z['ZYYR'].values
        Y22 = Z['ZYYI'].values
        
    else:
        freq = Z[:,0]
        nf = len(freq)
        
        ZXXR   = X11 = Z[:,1]
        ZXXI   = Y11 = Z[:,2]
        ZXYR   = X12 = Z[:,4]
        ZXYI   = Y12 = Z[:,5]
        ZYXR   = X21 = Z[:,7]
        ZYXI   = Y21 = Z[:,8]
        ZYYR   = X22 = Z[:,10]
        ZYYI   = Y22 = Z[:,11]


    #loop over each frequency 
    #for f in range(0, nf):
    det_x = (X11*X22 - X21*X12)
    ph11 = (X22*Y11 - X12*Y21)/det_x
    ph12 = (X22*Y12 - X12*Y22)/det_x 
    ph21 = (X11*Y21 - X21*Y11)/det_x
    ph22 = (X11*Y22 - X21*Y12)/det_x
    
    # Bibby et al. 2005 GJI
    pi1 = 0.5 * np.sqrt((ph11-ph22)**2 + (ph12+ph21)**2)
    pi2 = 0.5 * np.sqrt((ph11+ph22)**2 + (ph12-ph21)**2)
    
    phmax = np.degrees(np.arctan(pi2+pi1))
    phmin = np.degrees(np.arctan(pi2-pi1))
    
    ellip = (phmax-phmin) / (phmax+phmin)
    alpha = np.degrees(0.5 * np.arctan2((ph12+ph21) , (ph11-ph22)))
    beta  = np.degrees(0.5 * np.arctan2((ph12-ph21) , (ph11+ph22)))
    azimuth = alpha - beta
    
    ph_tens = np.zeros((nf, 5))
    ph_params = np.zeros((nf, 7))	
    ph_tens[:,0] = ph_params[:,0] = freq
    
    ph_tens[:,1] = ph11
    ph_tens[:,2] = ph12
    ph_tens[:,3] = ph21
    ph_tens[:,4] = ph22
    
    ph_params[:,1] = phmax
    ph_params[:,2] = phmin
    ph_params[:,3] = alpha
    ph_params[:,4] = beta
    ph_params[:,5] = ellip
    ph_params[:,6] = azimuth
    
    return(ph_tens, ph_params)
    



def phaseTensorErr(Z,MC_realizations = 10000):
    
    """
    Function that calculates phase tensor parameters errors
    It uses a Monte Carlo simulaiton 
    Following Booker (2014 Survey in Geophysics)
    
    
    Input:
    - data (pandas DF format for Z)
    - optional: number of MC realizations to perform (default =10000, should be enough)
    
    Output:
    - phase tensor matrix std. dev.
    - phase tensor parameters std. dev. (phmax,phmin,alpha,beta,ellip,azimuth)
    
    """   
        
    from scipy.stats import norm
    import sys

    
    freq = Z['FREQ'].values
    nf = len(freq)
    
    Z_MC = np.zeros((MC_realizations,12, len(freq)))
    #PT_MC = np.zeros((MC_realizations,2, len(freq)))
    ph_tens_MC = np.zeros((MC_realizations,4,len(freq)))
    ph_params_MC  = np.zeros((MC_realizations,6,len(freq)))
    
    ph_tens_std = np.zeros((len(freq),4))
    ph_params_std  = np.zeros((len(freq),6))
    #For each frequency:
    
    for f in range(nf):
    
        Z_MC[:,0,f] = freq[f]
        
        Z_MC[:,1,f] = np.random.normal(Z['ZXXR'][f], (Z['ZXX.VAR'][f])**0.5, MC_realizations)
        Z_MC[:,2,f] = np.random.normal(Z['ZXXI'][f], (Z['ZXX.VAR'][f])**0.5, MC_realizations)
        Z_MC[:,4,f] = np.random.normal(Z['ZXYR'][f], (Z['ZXY.VAR'][f])**0.5, MC_realizations)
        Z_MC[:,5,f] = np.random.normal(Z['ZXYI'][f], (Z['ZXY.VAR'][f])**0.5, MC_realizations)
        Z_MC[:,7,f] = np.random.normal(Z['ZYXR'][f], (Z['ZYX.VAR'][f])**0.5, MC_realizations)
        Z_MC[:,8,f] = np.random.normal(Z['ZYXI'][f], (Z['ZYX.VAR'][f])**0.5, MC_realizations)
        Z_MC[:,10,f] = np.random.normal(Z['ZYYR'][f], (Z['ZYY.VAR'][f])**0.5, MC_realizations)
        Z_MC[:,11,f] = np.random.normal(Z['ZYYI'][f], (Z['ZYY.VAR'][f])**0.5, MC_realizations)
    
    
        ph_tens, ph_params = phaseTensor(Z_MC[:,:,f])

        ph_tens_MC[:,:,f] = ph_tens[:,1:]
        ph_params_MC[:,:,f] = ph_params[:,1:]

        mean, ph_tens_std[f,0] = norm.fit(ph_tens_MC[:,0,f])
        mean, ph_tens_std[f,1] = norm.fit(ph_tens_MC[:,1,f])
        mean, ph_tens_std[f,2] = norm.fit(ph_tens_MC[:,2,f])
        mean, ph_tens_std[f,3] = norm.fit(ph_tens_MC[:,3,f])
        
        mean, ph_params_std[f,0] = norm.fit(ph_params_MC[:,0,f])
        mean, ph_params_std[f,1] = norm.fit(ph_params_MC[:,1,f])
        mean, ph_params_std[f,2] = norm.fit(ph_params_MC[:,2,f])
        mean, ph_params_std[f,3] = norm.fit(ph_params_MC[:,3,f])
        mean, ph_params_std[f,4] = norm.fit(ph_params_MC[:,4,f])
        mean, ph_params_std[f,5] = norm.fit(ph_params_MC[:,5,f])

    return(ph_tens_std, ph_params_std)



def medianFilter(freq, param, freq_sp):

    """
    Function that filters out outliers using a median filter ()
    """   
    
    param_filt = np.zeros(len(freq))
    # freq_sp = 0.2 # 1/2 interval value (log10) where the points can be considered for calculating the median
    for f1, f1_val in enumerate(freq):
        paramMed = []
        if f1 == 0:         # consider 0 beyond the beginning extremity (lim->0 = 0)
            paramMed = [0]
        for f2, f2_val in enumerate(freq):
            if np.log10(f1_val)-freq_sp < np.log10(f2_val) < np.log10(f1_val)+freq_sp:
                paramMed.append(param[f2])
        param_filt[f1] = np.median(np.array(paramMed))
    return param_filt




def z2rhophy(freq,Zr,Zi,dZ):
    
    import numpy as np

    FREQ = freq
    ZR = Zr
    ZI = Zi
    Z_VAR = dZ
    
    # calcul of apparent resistivity and phases
    rho = ((ZR**2+ZI**2)*0.2/(FREQ))
    phy = np.degrees(np.arctan2(ZI,ZR))
    
    # calcul of errors
    drho = ((ZR**2+ZI**2)**0.5)*(np.sqrt(Z_VAR))*0.4/(FREQ) 
    dphy= np.degrees(np.arcsin(np.sqrt(Z_VAR)/((ZR**2+ZI**2)**0.5)))

    return (rho, phy, drho, dphy)




def Zinv(Z):
    
    """
    Compute the invariant of Z: mean of the TE and TM modes
    """
    
    ZinvR = (Z['ZXYR'].values - Z['ZYXR'].values)/2
    ZinvI = (Z['ZXYI'].values - Z['ZYXI'].values)/2
    
    return ZinvR ,ZinvI 




def Zdet(Z):
    
    """
    Compute the determinant of Z
    
    Zdet = Zxx*Zyy - Zxy*Zyx = Zdet_1 - Zdet_2
    """
    nF = len(Z)
    # Calculate determinant 
    mat = np.zeros((2,2,nF), complex)
    mat[0,0,:] = Z['ZXXR'].values + (Z['ZXXI'].values * 1j)
    mat[0,1,:] = Z['ZXYR'].values + (Z['ZXYI'].values * 1j)
    mat[1,0,:] = Z['ZYXR'].values + (Z['ZYXI'].values * 1j)
    mat[1,1,:] = Z['ZYYR'].values + (Z['ZYYI'].values * 1j)
    
    ZdetR = np.zeros((nF))
    ZdetI = np.zeros((nF))
    sd = np.zeros((nF))
    lnSd = np.zeros((nF))
    
    for freq_det in range(len(Z)):
        det = np.linalg.det(mat[:,:,freq_det])**0.5
        ZdetR[freq_det] = det.real
        ZdetI[freq_det] = det.imag
    
    
    # Calculate determinant std dev.
    # (log transform of) x + dx --> log(x) + (log(1+dx/x) - log(1-dx/x))/(2*sqrt(2))
    cent = (ZdetR**2 + ZdetI**2)**0.5
    sd = (Z['ZXX.VAR'] + Z['ZXY.VAR'] + Z['ZYX.VAR'] + Z['ZYY.VAR'])**0.5
    lnSd = (np.log(1+sd/cent) - np.log(1-sd/cent))/(2*np.sqrt(2))
    
    return ZdetR ,ZdetI, sd, lnSd



def keepMax(x):
    # For the attributes we use the highest value encountered 
    # as a minimum value (see paper)
    for i in range(1,len(x)):
        if abs(x[i])>abs(x[i-1]):
            x[i] = abs(x[i])
        else:
            x[i] = abs(x[i-1])
    return x




def getData(edi_file_path, medfiltsp):
    
    # read edi file
    Z_orig, site_id, coord = readEDI(edi_file_path)
    
    # Conversion to appropriate units: we use ohms
    #   1 ohm = 10000(4*pi) [mV/km/nT]
    C = 10000/(4*np.pi)
    Z = Z_orig.copy(deep=True)
    Z.loc[:,['ZXXR','ZXXI','ZXYR','ZXYI','ZYXR','ZYXI','ZYYR','ZYYI']] = Z_orig.loc[:,['ZXXR','ZXXI','ZXYR','ZXYI','ZYXR','ZYXI','ZYYR','ZYYI']] / C
    Z.loc[:,['ZXX.VAR','ZXY.VAR','ZYX.VAR','ZYY.VAR']] = Z_orig.loc[:,['ZXX.VAR','ZXY.VAR','ZYX.VAR','ZYY.VAR']] / (C**2)

    freq = Z['FREQ']
    nF = len(freq)
    
    # calculate phase tensor parameters            
    ph_tens, ph_params = phaseTensor(Z)
    ph_tensERR, ph_paramsERR = phaseTensorErr(Z)
    beta = ph_params[:,4]
    ellip = ph_params[:,5]
    beta_err = ph_paramsERR[:,3]
    ellip_err = ph_paramsERR[:,4]

    # Calculate determinant of Z and error 
    ZdetR ,ZdetI, ZdetSd, ZdetLnSd = Zdet(Z)
    
    # Calculate differences between Zxy and Zyx
    difPol = abs(np.log(abs(Z['ZXYR']))-np.log(abs(Z['ZYXR']))) + abs(np.log(abs(Z['ZXYI']))-np.log(abs(Z['ZYXI'])))

    #Filter out outliers using a median filter
    ellip_filt = medianFilter(freq, ellip, medfiltsp)
    beta_filt = medianFilter(freq, beta, medfiltsp)
    difPol_filt = medianFilter(freq, difPol, medfiltsp)

    # For the attributes we use the highest value encountered as a minimum value (see paper)
    ellipM = keepMax(ellip_filt)
    betaM = keepMax(beta_filt)
    difPolM = keepMax(difPol_filt)

    # Create dataframes for inversion (data_inv) and general (data_Z)
    df_inv = pd.DataFrame(columns=['freq','ZdetRLn','ZdetILn','ZdetLnSd','elipM','betaM','difPolM'])
    df_inv['freq'] = freq
    df_inv['ZdetRLn'] = np.log(ZdetR)
    df_inv['ZdetILn'] = np.log(ZdetI)
    df_inv['ZdetLnSd'] = ZdetLnSd
    df_inv['elipM'] = ellipM
    df_inv['betaM'] = betaM
    df_inv['difPolM'] = difPolM
    
    dat = [Z_orig, ph_params, ph_paramsERR]
    
    return site_id, df_inv, dat 











def fwd1D(model, f, layers = 'depth'):

    
    """
    Compute the MT response of a 1D conductivity model
    
    Input:
    - model (nLayers x 2 array) 
        --> resistivity ohm.m
        --> depth in meters (last depth assumed to be the bottom of the model)
    - frequencies array at which the solution is computed (in Hz)
    
    Output:
    - impedance Z
    - apparent resistivity
    - phase
    
    """
    
    import numpy as np
    import math
   
    mu = 4*math.pi*1E-7; #Magnetic Permeability (H/m)
    
    res = model[:,0] 
    
	
    if layers == 'depth':
        depth = model[:,1]  #depth
        # Calculate thickness of each layers
        th = np.zeros((len(depth)))
        th[0] = depth[0]
        th[-1] = depth[-1]
        for i in range(1,len(depth[:-1])):
            th[i] = depth[i]-depth[i-1]
        #th[-2] =  th[-1] / 2
    else:
        th = model[:,0]  #thicknesses
	
    nlayers = len(res)
    
    con=1/res
    nf=len(f)
    
    ares=[]
    ares = np.array(ares, dtype = np.float32)
    phy=[]
    phy = np.array(phy, dtype = np.float32)
    
    #Define arrays
    w = 2*np.pi*f
    T = 1 / f
    
    Cm = np.empty(nlayers, complex)
    Z = np.empty(nf, complex)
    
    #Récurrence de Wait(1954)
    
    #Loop over each frequency
    for i in range(0, nf):    
        #1  Calculate bottom half space impedance
        gm = np.sqrt((w[i] * mu * (1.0/res[-1]))*1j);
        Cm[-1] = 1 / gm
        
        for k in range(nlayers-2, -1, -1):
            #2  Calculate impedance of layer k
            gm = np.sqrt((w[i] * mu * (1.0/res[k]))*1j);
            rm = (1 - gm*Cm[k+1]) / (1 + gm*Cm[k+1])
            Cm[k] = (1 -rm*np.exp(-2*gm*th[k])) / (gm*(1 +rm*np.exp(-2*gm*th[k])))
        Z[i] = Cm[0] *1j * w[i] * mu     #Cm[0] last Cm to be calculated: the one in surface
        # Units of Z in international standards units E/B(V/m/A/H)(EDI) from E/H(mV/km/nT) == factor mu
    
        # Step 3. Compute apparent resistivity and phase
        apparentResistivity = (abs(Z[i]) * abs(Z[i]))/(mu*w[i])
        phase = math.atan2(Z[i].imag, Z[i].real)
    
        ares=np.append(ares, apparentResistivity)
        phy=np.append(phy,phase)
    phy_deg=phy*180/math.pi #from radian to degrees
    
    return f, ares, phy_deg, Z




def niblettBostick_depthTransform(rho, phy, T):
    
    """
    "Immediate transformation of apparent resistivity and phase data and 
    presentation of an approximate resistivity and depth data. Specifically, 
    this transformation is based on the simple asymptotic expressions 
    introduced by Bostick (1977)"   
    
    --> GOLDBERG, S. and ROTSTEIN, Y. 1982, A Simple Form of Presentation of 
    Magnetotelluric Data Using the Bostick Transform,
    Geophysical Prospecting 30,211-216.

    Input:
    - apparent resistivity (Ohm meters)
    - phase angle (degrees)
    - period (seconds)
    
    Output:
    - resistivity estimate (Ohm meters)
    - depth (meters) 

    """
    
    mu0 = 4*np.pi*10**-7

    rho_nb = rho * (np.pi/(2*np.deg2rad(phy%90)) - 1)
    depth_nb = np.sqrt(rho*T/(2*np.pi*mu0))
    
    return rho_nb, depth_nb

    
