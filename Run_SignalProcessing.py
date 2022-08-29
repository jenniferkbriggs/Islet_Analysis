
# %% 
# This code is made to extract different frequencies of the signal and save them as individual csv files.
# Jennifer Briggs 2022

## ---- Options for you to change -------
# set to true if you'd like to see and save figures. set to false if you don't need figures
fig_on = True

# if you want to predefine a savepath. If not, comment out this line by putting at # in front!
global savepath
savepath = '/Users/jkbriggs/Documents/GitHub/Functional_and_Structural_Networks/Examples/'

lowcut = 0.1
highcut =  0.4


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

freqs = rfft(np.array(ca.iloc[:,3]))
plt.plot(freqs[0:1000])
# %% 

b, a = scipy.signal.butter(3, [lowcut, highcut], 'band', fs = fs)
filteredBandPass = scipy.signal.lfilter(b, a, ca)

plt.plot(filteredBandPass[:,1])
# %%
b, a = scipy.signal.butter(3, 0.000049, 'lowpass', fs = 1/fs)
filteredLowPass = scipy.signal.filtfilt(b, a, ca)

# b, a = scipy.signal.butter(3, 0.9, 'highpass', fs = fs)
# filteredHighPass = scipy.signal.filtfilt(b, a, ca)


plt.plot(time, ca.iloc[:,3])
plt.plot(time, filteredLowPass[:,3])
#plt.xlim([1,20])

#plt.plot(filteredHighPass[:,49])

# %%
windowSize = 400
window = np.hanning(windowSize)
window = window / window.sum()

# filter the data using convolution
filtered = np.convolve(window, ca.iloc[:,3], mode='valid')
plt.plot(filtered)
# %%

