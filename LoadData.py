#LoadData.py

# This is a function file which loads data. Here you can preset the path to give quick access to data or you can call a gui
# Jennifer Briggs - 2022. 

# Optional argument: path to timeseries file
# %%import packages
from Anne_code_spiketime_detection import OpenPrefined
import easygui
import pickle
import pandas as pd
import numpy as np

def LoadData(path = 0, USE_CONFIGURED_ISLETS = 'False', TEST_FILE):
    # if path is not passed in
    if path == 0:
        path = easygui.fileopenbox('Select Time signal file')

    #load file based on file type
    if path[-3:int(len(path))] == 'csv':
        ca = pd.read_csv(path, index_col=False)

    elif path[-3:int(len(path))] == 'mcr': #Anne's electrode data
        sd, td = OpenPrefined(path)
        #Electrode data is stored by x,y values. We will index them first along the x, then y. 
        # For example, [1,1] = 0, [1, 2] = 1, [2, 1] = len(y)
        numcell = (td.max_x - td.min_x+1)*(td.max_y - td.min_y+1)
        time = td.getElectrode(td.min_x, td.min_y)
        time = len(time.values)
        dat = np.empty([time, numcell])

        print('Reshaping Values')
        loc = np.empty([2,numcell])
        dat = dict()
        i = 0

        if USE_CONFIGURED_ISLETS:
            islet_configuration = read_islet_configuration(TEST_FILE)
            for islet in islet_configuration:
                for e in islet:
                    electrode = td.getElectrode(e[0],e[1])
                    #dat[:,i] = list(electrode.values)
                    dat.update({str(x) + ',' + str(y): list(electrode.values)})
                    loc[:,i] = [x,y]
                    i = i+1        
        else:
            for y in range(td.min_y, td.max_y+1):
                for x in range(td.min_x, td.max_x+1):
                    electrode = td.getElectrode(x,y)
                    #dat[:,i] = list(electrode.values)
                    dat.update({str(x) + ',' + str(y): list(electrode.values)})
                    loc[:,i] = [x,y]
                    i = i+1
        

        fs = td.tickrate 

        timeall = np.arange(0,int(time/fs),1/fs)
        ca = pd.DataFrame(dat)
        ca['Time'] = timeall


    print('Loaded Data Correctly')
    return ca 
# %%


# %%
