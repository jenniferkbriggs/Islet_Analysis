
# %% 
# This code is made to extract different frequencies of the signal and save them as individual csv files.
# Jennifer Briggs 2022

## ---- Options for you to change -------
# set to true if you'd like to see and save figures. set to false if you don't need figures
fig_on = True

# if you want to predefine a savepath. If not, comment out this line by putting at # in front!
global savepath
savepath = '/Users/jkbriggs/Documents/GitHub/Functional_and_Structural_Networks/Examples/'



# %% Import packages
from SignalProcessing import *
import pandas as pd
import easygui
import numpy as np


# %%  Load calcium file
try: # you can directly add the path to your data here
    ca = pd.read_csv('/Users/jkbriggs/OneDrive - The University of Colorado Denver/Anschutz/Islet/TempData/Erli_calcium.csv')
except: #if there is no path, it will ask you to select the folder
    path = easygui.fileopenbox('Select Time signal file')
    ca = pd.read_csv(path)

try: #if time is in the first axis, we save it and remove
    time = ca.Time
    ca = ca.drop('Time', axis=1)
    tdiff = list(time.diff())
    fs = 1/np.nanmean(tdiff) #note that the value of this is not Hz but 1/ what ever time interval used (usually ms)
except:
    print('No time avaliable')
    fs = input('What is the frequency of recording?')
# %%
