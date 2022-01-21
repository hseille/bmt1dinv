'''
    File name: DDE.py
    Authors: Hoël Seillé / Gerhard Visser
    Date created: 01/06/2020
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
# Define here the files and parameters to use to create the input data
# =============================================================================

# specify project folder name
project = 'example'

# error floor to use, NOT in percent! 
# -1 indicates that no error floors are used, just dimensionality errors
EF = -1

# it defines the lenght of the window, that will be used to perform the
#    median filter, in log scale, along the frequency axis
# medfiltsp = 0.5 => median filter along half a frequency decade
medfiltsp = 0.5

# options to plot the data and create the csv files needed for the inversion
plotMTdata = True
saveCSVfiles = True

# option to consider possible static shift: if StSh is set to True, an error
#   equivalent to the split between XY and YX will be assigned. THis option
#   will cause small but potentially relevant signal to be ignored.
#   default = False
StSh = False

# Define the path to the error model files used, it shouldn't be changed 
#   unless a new tree is defined
atts_file = '../projects/%s/tree/atts.txt'%(project)
DDMfile  =  '../projects/%s/tree/tree.txt'%(project)

edi_path  = '../projects/%s/edi'%(project)
csv_path  = '../projects/%s/csv'%(project)
# =============================================================================
# 
# =============================================================================










# Load modules
import os
import numpy as np
import warnings
warnings.filterwarnings('ignore')

import tree as tree
import MT as MT
import plots as plots


print('\n#############################################')
print('## 1-D Magnetotellurics Bayesian Inversion ##')
print('##  Dimensionality Discrepancy Estimation  ##')
print('##       Hoël Seillé / Gerhard Visser      ##')
print('##          CSIRO DEI FSP   2020           ##')
print('#############################################\n')





#read and store the tree attribute list 
atts_dic = tree.readAttsFile(atts_file)

#read and store the decision tree 
tr = tree.readTreeFile(DDMfile, atts_dic)

print('EDI files folder: ',edi_path)
# loop over the data files (EDI files)
for root, dirs, files in os.walk(edi_path):
    for file in files:
        if file.endswith('.edi'):
            edi_file = file
            print('\nMT station %s... '%edi_file)
            edi_file_path = r'%s/%s'%(root, edi_file)
            print('    Loading EDI file... ')

            # Create Data file
            site_id, data_1D, dataMT = MT.getData(edi_file_path, medfiltsp)
            
            # Assess data_1D
            if data_1D.isnull().values.any():
                bad_f = np.where(np.isnan(data_1D))[0]
                for fff in range(len(bad_f)):
                    print('    .. invalid value found for freq number %d (f = %5.5f Hz). Removed.'%(bad_f[fff], data_1D['freq'][bad_f].values[fff]))
                data_1D = data_1D.drop(bad_f)
                data_1D = data_1D.reset_index(drop=True)
                dataMT[0] = (dataMT[0]).drop(bad_f)
                dataMT[0] = (dataMT[0]).reset_index(drop=True)
                dataMT[1] = np.delete(dataMT[1],bad_f,0)
                dataMT[2] = np.delete(dataMT[2],bad_f,0)

            # Get covariance matrix
            C, ss = tree.getC(data_1D, tr)

            # load and plot MT data
            if plotMTdata: 
                import matplotlib.pyplot as plt
                plots.plot_edi(site_id,data_1D, ss, dataMT, medfiltsp, 
                               plot_rhoPhy=True, 
                               plot_antidiag=True, 
                               plot_diag=True,
                               phase_90=True, 
                               plot_size = 'small')
                plt.savefig('%s/%s.png'%(edi_path,site_id),dpi=600, 
                            bbox_inches="tight")
                plt.close('all')
                print('    Data plot saved');
                
            #  save CSV files for the inversion
            if saveCSVfiles:
                if EF < 0:
                    df,ss = tree.exportForTransD(site_id,data_1D, 
                                                 tr, errorfloor=-1)
                    df.to_csv(r'%s/%s.csv'%(csv_path,site_id),
                              sep=',', index=False)

                else:
                    df,ss = tree.exportForTransD(site_id,data_1D, 
                                                 tr, errorfloor=EF)
                    df.to_csv('%s/%sEF%d.csv'%(csv_path,site_id,100*EF),
                                                  sep=',', index=False)

                print('    Input .csv file saved');

