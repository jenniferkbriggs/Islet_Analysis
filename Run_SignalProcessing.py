
# %% 
# This code is made to extract different frequencies of the signal and save them as individual csv files.
# Jennifer Briggs 2022

## ---- Options for you to change -------
# set to true if you'd like to see and save figures. set to false if you don't need figures
fig_on = True

# Define window size
windowSize = 4000

# if you want to predefine a savepath. If not, comment out this line by putting at # in front!
global savepath

lowcut = 1/10000
highcut =  100


# %% Import packages
from SignalProcessing import *
import pandas as pd
import easygui
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import freqz
from LoadData import *
from scipy.fft import rfft,rfftfreq, irfft
import scipy


# %%  Load calcium file
path = '/Users/jkbriggs/Dropbox/CMOS data/210720_3985_G10.cmcr'
savepath = '/Users/jkbriggs/Documents/GitHub/Islet_Analysis/Examples/SignalProcessing/'
ca = LoadData(path)

# %%  Load calcium file
try: #if time is in the first axis, we save it and remove
    time = ca.Time
    ca = ca.drop('Time', axis=1)
    tdiff = list(time.diff())
    fs = 1/np.nanmean(tdiff) #note that the value of this is not Hz but 1/ what ever time interval used (usually ms)
except:
    print('No time avaliable')
    fs = int(input('What is the frequency of recording?'))

# %%
window = np.hanning(windowSize)
window = window / window.sum()

# filter the data using convolution
slowar = np.apply_along_axis(lambda m: np.convolve(window, m, mode='same'), axis=0, arr=ca)
fast = ca - slowar

# %%
slow = pd.DataFrame(slowar)
slow["Time"] = time
fast["Time"] = time

#save
if 'savepath' not in locals():
    savepath = easygui.diropenbox('Select Folder to Save Data In')
    savepath = savepath + savepath[0] #adds slash
try: 
    slow.to_csv(savepath + 'Slow.csv', index=False)
    fast.to_csv(savepath + 'Fast.csv', index=False)
except:
    print('Save path is not working')
    savepath = easygui.diropenbox('Select Folder to Save Data In')
    savepath = savepath + savepath[0] #adds slash
    slow.to_csv(savepath + 'Slow.csv', index=False)
    fast.to_csv(savepath + 'Fast.csv', index=False)

if fig_on == True: #save and close figures
    plt.plot(time, slowar[:,3])
    plt.title('Low Frequency')
    plt.xlabel('Time')
    plt.ylabel('Voltage')
    plt.savefig(savepath + 'LowFreq.png')
    plt.clf 

    plt.plot(time, fast.iloc[:,3])
    plt.title('High Frequency')
    plt.xlabel('Time')
    plt.ylabel('Voltage')
    plt.savefig(savepath + 'HighFreq.png')
    plt.clf #close

print('Done')
# %% Forier Filtering
# freqs = rfft(np.array(ca.iloc[:,3]))
# plt.plot(freqs[0:1000])
# # %% 

# b, a = scipy.signal.butter(3, [lowcut, highcut], 'band', fs = fs)
# filteredBandPass = scipy.signal.lfilter(b, a, ca)
# plt.plot(time, ca.iloc[:,4])
# plt.plot(time, filteredBandPass[:,4])


# %% 
#df = pd.DataFrame(filteredBandPass)
#plt.plot(filteredHighPass[:,49])

