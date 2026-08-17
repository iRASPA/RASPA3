import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pathlib
import os
from os import path

from scipy.signal import argrelextrema

def GetMinMax(Dataset, colname, p0):
    max_ind = argrelextrema(np.array(Dataset[colname]), np.greater)
    rmax    = np.array(Dataset[colname])[max_ind]
    min_ind = argrelextrema(np.array(Dataset[colname]), np.less)
    rmin    = np.array(Dataset[colname])[min_ind]
    gas_N   = np.min(max_ind)
    gas_lnpi= np.array(Dataset[colname])[gas_N]
    liq_N   = np.max(max_ind)
    liq_lnpi= np.array(Dataset[colname])[liq_N]
    sads    = (min_ind > gas_N) & (min_ind < liq_N)

    print(liq_N, rmax)

    # Get the minimum among the saddles #
    min_saddles_idx = list(rmin).index(np.min(rmin))
    #min_saddles_N = min_ind[0][min_saddles_idx]
    #print(min_saddles_idx, min_ind[0][min_saddles_idx])
    #sad_N   = np.min(np.array(min_ind)[sads])
    sad_N   = min_ind[0][min_saddles_idx]
    sad_lnpi= np.array(Dataset[colname])[sad_N]
    
    
    return [p0, sad_N, sad_lnpi, gas_N, gas_lnpi, liq_N, liq_lnpi]

def Combine_TMMC_lnpi(file, file_dir, start, sims, p0):
    print(file_dir + str(start) + '/' + file)
    TMMCFile = file_dir + str(start) + '/' + file
    if(path.exists(TMMCFile)):
        InitData = pd.read_table(TMMCFile, sep='\s+', 
                                 skiprows = 12)#, skipfooter = end)
        MinMaxData = pd.DataFrame(columns=['Fugacity', 'Saddle-N', 'Saddle-lnpi',
                                           'gas-N', 'gas-lnpi', 
                                           'liq-N', 'liq-lnpi'])
        for sim in range(start + 1, sims):
            TMMCFile = file_dir + str(sim) + '/' + file
            if(path.exists(TMMCFile)):
                Dataset = pd.read_table(TMMCFile, sep='\s+', 
                                        skiprows = 12)#, skipfooter = end)
    
                diff_lnpi = Dataset['lnpi'] - Dataset['lnpi'][0]
                last_lnpi = InitData.iloc[-1]['lnpi']
                Dataset['lnpi'] = last_lnpi + diff_lnpi
                # Drop first row
                Dataset.drop(index=Dataset.index[0], axis=0, inplace=True)
                frames = [InitData, Dataset]
                InitData = pd.concat(frames)
        # Normalize
        max_lnpi = np.max(InitData['lnpi'])
        print("Max_lnpi: " + str(max_lnpi))
        InitData['cont-maxlnpi'] = InitData['lnpi'] - max_lnpi
        InitData['pi'] = np.exp(InitData['cont-maxlnpi'])
        sumpi = np.sum(InitData['pi'])
        print("SumPi: " + str(sumpi))
        InitData['normalized-pi'] = InitData['pi'] / sumpi
        InitData['normalized-lnpi'] = np.log(InitData['normalized-pi'])
        
        row = GetMinMax(InitData, 'normalized-lnpi', p0)
        MinMaxData.loc[len(MinMaxData)] = row
        lnpis = [row[2], row[4], row[6]]
        min_lnpi = np.min(lnpis)
        max_lnpi = np.max(lnpis)
        
        figpX = plt.figure()
        axpX = figpX.add_subplot(111)
        axpX.plot(InitData['N'], InitData['normalized-lnpi'], color='orange',
                  label='lnpi', linewidth = aaa)
        axpX.set_xlabel('Loading [Number of Molecule]', fontsize = fsize, fontweight='bold')
        axpX.set_ylabel("lnpi", fontsize = fsize, fontweight='bold')
        #axpX.set_ylim([init_lnpi,max_lnpi])
        axpX.set_ylim([min_lnpi - 10, max_lnpi + 10])
        plt.savefig(file_dir + '/' + 'combined_lnpi' + '.png', dpi=900, bbox_inches = 'tight')
        plt.clf()
        
        return InitData, MinMaxData
    
def analyse_TMMC_lnpi(file, file_dir):
    if(path.exists(file)):
        Dataset = pd.read_table(TMMCfile, sep='\s+', 
                                skiprows = 12)#, skipfooter = end)
        figpX = plt.figure()
        axpX = figpX.add_subplot(111)
        axpX.plot(Dataset['N'], Dataset['lnpi'], color='orange',
                  label='lnpi', linewidth = aaa)
        axpX.set_xlabel('Loading [Number of Molecule]', fontsize = fsize, fontweight='bold')
        axpX.set_ylabel("lnpi", fontsize = fsize, fontweight='bold')
        plt.savefig(file_dir + '/' + 'lnpi_Sim_' + str(sim) + '.png', dpi=900, bbox_inches = 'tight')
        plt.clf()

def Reweight_TMMC_lnpi(Dataset, MinMaxData, p0, p1):
    Dataset['lnpi_' + str(p1)] = Dataset['normalized-lnpi'] + Dataset['N'] * np.log(p1/p0)
    
    # Add normalization
    max_lnpi = np.max(Dataset['lnpi_' + str(p1)])
    cont_maxlnpi = Dataset['lnpi_' + str(p1)] - max_lnpi
    cont_maxpi = np.exp(cont_maxlnpi)
    sumpi = np.sum(cont_maxpi)

    normalized_pi = cont_maxpi / sumpi
    Dataset['lnpi_' + str(p1)] = np.log(normalized_pi)
    
    row = GetMinMax(Dataset, 'lnpi_' + str(p1), p1)
    MinMaxData.loc[len(MinMaxData)] = row
    
    lnpis = [row[2], row[4], row[6]]
    min_lnpi = np.min(lnpis)
    max_lnpi = np.max(lnpis)
    
    figpX = plt.figure()
    axpX = figpX.add_subplot(111)
    axpX.plot(Dataset['N'], Dataset['lnpi_' + str(p1)], color='orange',
              label='lnpi', linewidth = aaa)
    init_lnpi = Dataset.iloc[0]['lnpi_' + str(p1)]
    end_lnpi  = Dataset.iloc[-1]['lnpi_' + str(p1)]
    selected_lnpi = np.max([init_lnpi, end_lnpi])
    max_lnpi  = np.max(Dataset['lnpi_' + str(p1)])
    #axpX.set_ylim([init_lnpi,max_lnpi])
    axpX.set_ylim([min_lnpi - 10, max_lnpi + 10])
    axpX.set_xlabel('Loading [Number of Molecule]', fontsize = fsize, fontweight='bold')
    axpX.set_ylabel("lnpi", fontsize = fsize, fontweight='bold')
    plt.savefig(file_dir + '/' + 'lnpi_reweightedTo_' + str(p1) + '.png', dpi=900, bbox_inches = 'tight')
    plt.clf()

    return Dataset, MinMaxData

def GetMaxMin(Dataset, p1):
    return Dataset

fsize = 15
legendwidth = 4.0
aaa = 2
sim = 5

file_dir = os.getcwd() + '/'

TMMCfile = 'tmmc/tmmc_statistics.txt'
FinalData, MinMaxData = Combine_TMMC_lnpi(TMMCfile, file_dir, 0, 5, 36000)
# preslist = [35000,35500,36000,36500,37000,37500,38000,38500,39000,39500,40000,40500,41000,
#             41500,42000,42500,43000,43500,44000,44500,45000,45500,
#             46000,46500,47000,47500,48000,48500,49000,49500,50000,50500,51000,51500,52000]
preslist = np.linspace(42000 + 500, 50000, 16)

#preslist = [35000]
for pres in preslist:
    FinalData, MinMaxData = Reweight_TMMC_lnpi(FinalData, MinMaxData, 36000, pres)
FinalData.to_csv("combined_lnpi.csv")
MinMaxData.to_csv("combined_min_max-lnpi.csv")
# FinalData, MinMaxData = Reweight_TMMC_lnpi(FinalData, 36000, 45053.3)
# column = 'lnpi_' + '45053.3'
# max_ind = argrelextrema(np.array(FinalData[column]), np.greater)
# rmax    = np.array(FinalData[column])[max_ind]
# min_ind = argrelextrema(np.array(FinalData[column]), np.less)
# rmin    = np.array(FinalData[column])[min_ind]
