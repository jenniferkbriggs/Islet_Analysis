#region imports
import sys, importlib, os
import tkinter.filedialog
from builtins import enumerate

import McsPy.McsData
import McsPy.McsCMOS
from McsPy import ureg, Q_, McsCMOSMEA

import matplotlib
import seaborn as sns

import pickle
import scipy as sp
import scipy.signal
import numpy as np
from tkinter import *
import tkinter as tk
import xlsxwriter
from IPython.display import HTML
from ipywidgets import *
from matplotlib import pyplot as plt
from operator import itemgetter
#endregion

plt.style.use('ggplot')
#region Classes



def OpenPrefined(path):
    global td
    global sd
    global filename
    data = McsCMOSMEA.McsData(path)
    sd = data.Acquisition.Sensor_Data
    td = getTestData(sd)
    x_values = range(td.min_x,td.max_x+1)
    y_values = range(td.min_y,td.max_y+1)

    return sd, td




class Testdata:
    def __init__(self, electrodes, dictionary, tickrate, min_x, max_x, min_y, max_y):
        self.electrodes = electrodes
        self.dictionary = dictionary
        self.tickrate = tickrate        #tickrate in microseconds
        self.min_x = min_x
        self.max_x = max_x
        self.min_y = min_y
        self.max_y = max_y

    def getElectrode(self,x,y):
        return self.electrodes[self.dictionary[x][y]]

    def getElectrodeId(self,x,y):
        return self.dictionary[x][y]

    def getElectrodeNeighbours(self, electrode):
        neigh = []
        for i in [-1,0,1]:
            for j in [-1,0,1]:
                if self.getElectrodeId(electrode.x + i,electrode.y + j) is not None:
                    if i != 0 or j != 0:
                        neigh.append(self.getElectrode(electrode.x + i,electrode.y + j))
        return neigh

    def getMaximumPeakScore(self):
        ps = 0
        for el in self.electrodes:
            _ps = el.getMaximumPeakScore()
            if _ps > ps:
                ps = _ps
        return ps

    def getSaveData(self):
        el = []
        for e in self.electrodes:
            el.append(e.getSaveData())
        sd = Savedata(el, self.min_x, self.max_x, self.min_y, self.max_y)
        return sd

class Savedata:
    def __init__(self, electrodes, min_x, max_x, min_y, max_y):
        self.electrodes = electrodes
        self.min_x = min_x
        self.max_x = max_x
        self.min_y = min_y
        self.max_y = max_y

class Electrode:
    def __init__(self, x, y, conversion_factor, values=[], peaks=[]):
        self.x = x
        self.y = y
        self.values = values
        self.conversion_factor = conversion_factor
        self.peaks = peaks
        self.signatures = []
        self.matches = [[None for y in range(3)] for x in range(3)]

    def addMatches(self, matches, electrode):
        x = electrode.x - self.x # -1 if left 1 if right ->
        y = electrode.y - self.y # -1 if above 1 if below v
        if x >= -1 and x <= 1 and y >= 1 and y <= 1:
            self.matches[1-x][1-y] = matches

    def getMatches(self, electrode):
        x = electrode.x - self.x  # -1 if left 1 if right ->
        y = electrode.y - self.y  # -1 if above 1 if below v
        if x >= -1 and x <= 1 and y >= 1 and y <= 1:
            return self.matches[1 - x][1 - y]
        else:
            return None

    def getSaveData(self):
        return Electrode(self.x, self.y, self.conversion_factor, peaks=self.peaks)

    def getPeak(self, timestamp):
        for peak in self.peaks:
            if peak.timestamp == timestamp:
                return peak

    def getAmountBigPeaks(self, threshold):
        amount = 0
        for peak in self.peaks:
            if peak.peakScore > threshold:
                amount = amount + 1
            if (peak.peakScore * (-1)) > threshold:
                amount = amount + 1
        return amount

    def getMaximumPeakScore(self):
        ps = 0
        for peak in self.peaks:
            if peak.peakScore > ps:
                ps = peak.peakScore
        return ps

    def getMinimumPeakScore(self):
        ps = 0
        for peak in self.peaks:
            if peak.peakScore < ps:
                ps = peak.peakScore
        return ps



    # def getShowPeakPositions(self, show):
    #     positions = []
    #     if(show == 'all'):
    #         for peak in self.peaks:
    #             positions.append(peak.position)
    #     elif(show == 'maximums'):
    #         for peak in self.peaks:
    #             if(peak.isMaximum == 1):
    #                 positions.append(peak.position)
    #     elif (show == 'minimums'):
    #         for peak in self.peaks:
    #             if (peak.isMaximum == 0):
    #                 positions.append(peak.position)
    #
    #     return positions

    def getPeakPositions(self, threshold = 0):
        positions = []
        for peak in self.peaks:
            if peak.peakScore > threshold:
                positions.append(peak.position)
            if (peak.peakScore * (-1)) > threshold:
                positions.append(peak.position)

        return positions

    def getPeakTimestamps(self):
        timestamps = []
        for peak in self.peaks:
            timestamps.append(peak.timestamp)
        return timestamps

    def getPeakDeviations(self):
        deviations = []
        for peak in self.peaks:
            deviations.append(peak.deviation)
        return deviations

    def getPeakValues(self):
        values = []
        for peak in self.peaks:
            values.append(peak.value)
        return values

class Peak:
    def __init__(self, value, position, peakScore, timestamp):
        self.value = value
        self.position = position
        self.peakScore = peakScore
        self.timestamp = timestamp
        #self.isMaximum = isMaximum # 1 if Maximum, 0 if Minimum
        self.normalizedPeakScore = 0
        self.normalizedPosition = 0

    def distanceTo(self, otherPeak):
        vDistance = abs(self.position - otherPeak.position)
        hDistance = abs(self.peakScore - otherPeak.peakScore)
        return hDistance+vDistance

    def getWeight(self):
        return 1/3

#endregion

#region Global variables

global sd
sd = None
global td
td = None
global time_conversion_factor
time_conversion_factor = 0.000001
global filename
filename = ''

#endregion

#region Global functions

def getTestData(sd):
    global td
    electrodes = []
    # 65x65 with offset 1
    dictionary = [[None for y in range(67)] for x in range(67)]
    sm = sd.SensorMeta
    min_y = sm['Region.Top'][0]
    min_x = sm['Region.Left'][0]
    max_y = sm['Region.Bottom'][0]
    max_x = sm['Region.Right'][0]
    if (len(sm) == 1):
        #Alle Elektroden
        sdata = sd.SensorData_1_1
        counter = 0
        #Format [time][x][y]
        conversionFactors = sd.SensorMeta['Conversion Factors'][0]
        conversionFactors = conversionFactors.split(b' ')
        #y and x are swapped in the raw data in this case for no apparent reason
        for x in range(sdata.shape[1]):
            for y in range(sdata.shape[2]):
                #print('x: ', x+min_x , ' y: ', y+min_y)
                conversionFactor = int(conversionFactors[counter]) / 1000
                electrodes.append(Electrode(x+min_x,y+min_y,conversionFactor, sdata[:,x,y]))
                dictionary[x+min_x][y+min_y] = counter
                counter+=1
    else:
        #Ausgeählte Elektroden
        #atts = [a for a in dir(sd) if a.startswith('SensorData_') and not callable(getattr(sd, a))]

        counter = 0
        min_x = np.amin(sm['Region.Left'][:])
        max_x = np.amax(sm['Region.Left'][:])
        min_y = np.amin(sm['Region.Top'][:])
        max_y = np.amax(sm['Region.Top'][:])
        for i in range(len(sm)):
            att_name = 'SensorData_' + str(i + 1) + '_1'
            values = getattr(sd, att_name)
            x = sm['Region.Left'][i]
            y = sm['Region.Top'][i]
            conversionFactor = int(sm['Conversion Factors'][i]) / 1000
            #print('x: ', x, ' y: ', y)
            electrodes.append(Electrode(x, y, conversionFactor, values[:,0,0]))
            dictionary[x][y] = counter
            counter +=1
            #print(values[len(values) - 1])
            # print('number: ', i, ' Top: ', sm['Region.Top'][i], ' Left: ', sm['Region.Left'][i], ' Bottom: ', sm['Region.Bottom'][i], ' Right: ', sm['Region.Right'][i])
    td = Testdata(electrodes, dictionary, sd.SensorMeta['Tick'][0], min_x, max_x, min_y, max_y)
    print("Loading successful. loaded " + str(len(td.electrodes)) + " electrodes with " + str(len(td.electrodes[0].values)) + " values each.")
    return td

#endregion

#region tkInter Button functions

def refreshDropdown(dropdown, dd_var, newList):
    dropdown['menu'].delete(0,'end')
    for item in newList:
        dropdown['menu'].add_command(label=str(item), command=tk._setit(dd_var, item))

def refreshPeaks(*args):
    x = int(x_dropdown_variable.get())
    y = int(y_dropdown_variable.get())
    global td
    electrode = td.getElectrode(x, y)
    peaks = electrode.getPeakTimestamps()

def onFilemenuOpen():
    global td
    global sd
    global filename
    loading_label.config(text="Loading...")
    myFile = tkinter.filedialog.askopenfile(filetypes=[("CMOS MEA Files", "*.cmcr")])
    data = McsCMOSMEA.McsData(myFile.name)
    filename = os.path.basename(myFile.name)
    sd = data.Acquisition.Sensor_Data
    td = getTestData(sd)
    x_values = range(td.min_x,td.max_x+1)
    refreshDropdown(x_dropdown,x_dropdown_variable,x_values)
    y_values = range(td.min_y,td.max_y+1)
    refreshDropdown(y_dropdown, y_dropdown_variable, y_values)
    loading_label.config(text="Loading successful. Loaded " + str(len(td.electrodes)) + " electrodes with " + str(len(td.electrodes[0].values)) + " values each.")
    file_label.config(text="File : " + filename)



def onCalculatePeaks():
    x = int(x_dropdown_variable.get())
    y = int(y_dropdown_variable.get())
    method = method_dropdown_variable.get()
    global td
    tickrate = td.tickrate
    electrode = td.getElectrode(x, y)
    window_size_average = int(window_size_average_text.get("1.0", tk.END).strip('\n'))
    window_size_argrel = int(window_size_argrel_text.get("1.0", tk.END).strip('\n'))
    tolerance = float(std_text.get("1.0", tk.END).strip('\n'))

    #calculateMaximums(electrode,window_size_average,tolerance,window_size_argrel)


    electrode.peaks = []
    if(method == "Average Distance"):
        ps = calculatePeakScorePeaksAD(electrode.values, window_size_average, tolerance, window_size_argrel)
    if(method == "Maximum Distance"):
        ps = calculatePeakScorePeaksMD(electrode.values, window_size_average, tolerance, window_size_argrel)
    for p in ps:
        electrode.peaks.append(Peak(electrode.values[p[0]], p[0], p[1], float(p[0])*float(tickrate)))

    # peaks = maxs
    # peaks.extend(mins)
    # peaks.sort()
    # for p in peaks:
    #     electrode.peaks.append(Peak(electrode.values[p],p,0,float(p)*float(tickrate)))

    peaks = electrode.getPeakTimestamps()

def onCalculateAllPeaks():
    global td
    tickrate = td.tickrate
    # x = int(x_dropdown_variable.get())
    # y = int(y_dropdown_variable.get())
    method = method_dropdown_variable.get()
    window_size_average = int(window_size_average_text.get("1.0", tk.END).strip('\n'))
    window_size_argrel = int(window_size_argrel_text.get("1.0", tk.END).strip('\n'))
    tolerance = float(std_text.get("1.0", tk.END).strip('\n'))

    for electrode in td.electrodes:
        if (method == "Average Distance"):
            ps = calculatePeakScorePeaksAD(electrode.values, window_size_average, tolerance, window_size_argrel)
        if (method == "Maximum Distance"):
            ps = calculatePeakScorePeaksMD(electrode.values, window_size_average, tolerance, window_size_argrel)
        #calculateMaximums(electrode, window_size_average, tolerance, window_size_argrel)
        electrode.peaks=[]
        for p in ps:
            electrode.peaks.append(Peak(electrode.values[p[0]],p[0],p[1],float(p[0])*float(tickrate)))
    # electrode=td.getElectrode(x,y)
    # peaks = electrode.getPeakTimestamps()
    # refreshDropdown(peak_dropdown, peak_dropdown_variable, peaks)


def onShowElectrodePlot():
    global td
    x = int(x_dropdown_variable.get())
    y = int(y_dropdown_variable.get())
    threshold = float(threshold_big_text.get('1.0', 'end-1c'))

    electrode = td.getElectrode(x,y)
    tr = td.tickrate
    xs = range(len(electrode.values))
    xs = xs * tr
    xs = xs / 1000000
    fig, ax = plt.subplots()
    markers_on = electrode.getPeakPositions(threshold)
    ax.plot(xs, electrode.values, '-o', markevery=markers_on, color='tab:blue', mfc='r')
    ax.set_xlabel('time in s')
    ax.set_ylabel('μV')
    ax.set_title('X: ' + str(electrode.x) + ' Y: ' + str(electrode.y))
    plt.show()

#TODO:
def onShowPeakFrequencyHeatmap():
    global td
    # seconds = len(td.electrodes[0].values) * td.tickrate / 1000000
    threshold = float(threshold_big_text.get('1.0', 'end-1c'))
    frequencies = np.zeros((66,66))
    for electrode in td.electrodes:
        #peaks = len(electrode.peaks)
        peaks = electrode.getAmountBigPeaks(threshold)
        frequencies[electrode.y,electrode.x] = peaks # / seconds

    sns.heatmap(frequencies, xticklabels=True, yticklabels=True, annot=True, fmt='.5g')
    plt.xlim(0,66)
    plt.ylim(0,66)
    plt.gca().invert_yaxis()
    plt.show()


def onCreateExcelFile():
    global td

    global filename
    threshold = float(threshold_big_text.get('1.0', 'end-1c'))

    saveName = filename.split(".")
    myFile = tkinter.filedialog.asksaveasfile(mode='w', defaultextension=".xlsx", initialfile=saveName[0])

    workbook = xlsxwriter.Workbook(myFile.name)
    worksheet = workbook.add_worksheet()
    worksheet.write(0,0,'Y\X')
    for i in range(1,66):
        worksheet.write(0,i,i)
        worksheet.write(i,0,i)
    for el in td.electrodes:
        worksheet.write(el.y, el.x, el.getAmountBigPeaks(threshold))
    workbook.close()

def onSaveFile():
    #TODO: Ensure pkl Datatype
    global td
    global filename

    saveName = filename.split(".")
    myFile = tkinter.filedialog.asksaveasfile(mode='w', defaultextension=".pkl",initialfile=saveName[0])

    with open(myFile.name, 'wb') as outp:
        pickle.dump(td, outp, pickle.HIGHEST_PROTOCOL)

def onLoadFile():
    #TODO: prüfe ob richtige TestData
    global td
    global filename

    myFile = tkinter.filedialog.askopenfile(filetypes=[("pkl Files", "*.pkl")])
    filename=os.path.basename(myFile.name)

    with open(myFile.name, 'rb') as inp:
        saveData = pickle.load(inp)
        td = saveData
    x_values = range(td.min_x, td.max_x + 1)
    refreshDropdown(x_dropdown, x_dropdown_variable, x_values)
    y_values = range(td.min_y, td.max_y + 1)
    refreshDropdown(y_dropdown, y_dropdown_variable, y_values)
    loading_label.config(text="Loading successful. Loaded " + str(len(td.electrodes)) + " electrodes with " + str(
        len(td.electrodes[0].values)) + " values each.")
    file_label.config(text="File : " + os.path.basename(myFile.name))
    refreshPeaks()

#Calculate a Peak Score for each pre-filtered local extrema. Peak Score is calculated by average distance around extrema
def calculatePeakScorePeaksAD(values,k,h,w):
    peakScores = []
    peakScoresMax = []#
    peakPositionsMax = []
    peakScoresMin = []
    peakPositionsMin = []
    #filter_window_size = int(window_size_argrel_text.get("1.0", tk.END).strip('\n'))
    #h = float(std_text.get("1.0", tk.END).strip('\n'))
    #k = int(window_size_average_text.get("1.0", tk.END).strip('\n'))
    max_type = max_dropdown_variable.get()
    maximums = scipy.signal.argrelextrema(values, np.greater, 0, w)
    minimums = scipy.signal.argrelextrema(values, np.less, 0, w)

    if max_type == "All" or max_type == "Maximums":
        for point in maximums[0]:
            if point > k:
                if point < len(values) - k:
                    value = values[point]
                    signed_top = 0
                    signed_bot = 0
                    for i in range(k):
                        signed_bot += value - values[point-i]
                        signed_top += value - values[point+i]
                    ps = (signed_top / k + signed_bot / k ) / 2
                    if(ps > 0):
                        peakScoresMax.append(ps)
                        peakPositionsMax.append(point)

        std_max = np.std(peakScoresMax)
        mean_max = np.average(peakScoresMax)
        for i in range(len(peakScoresMax)):
            if(peakScoresMax[i] - mean_max > h * std_max):
                peakScores.append([peakPositionsMax[i], peakScoresMax[i]])

    if max_type == "All" or max_type == "Minimums":
        for point in minimums[0]:
            if point > k:
                if point < len(values) - k:
                    value = values[point]
                    signed_top = 0
                    signed_bot = 0
                    for i in range(k):
                        signed_bot += value - values[point-i]
                        signed_top += value - values[point+i]
                    ps = (signed_top / k + signed_bot / k ) / 2
                    if(ps < 0) :
                        peakScoresMin.append(ps)
                        peakPositionsMin.append(point)

        std_min = np.std(peakScoresMin)
        mean_min = np.average(peakScoresMin)
        for i in range(len(peakScoresMin)):
            if (peakScoresMin[i] - mean_min < h * std_min * -1):
                peakScores.append([peakPositionsMin[i], peakScoresMin[i]])


    peakScores.sort(key=itemgetter(0))

    return peakScores

#Calculate a Peak Score for each data point. Peak Score is calculated by average distance around extrema
def calculatePeakScorePeaks(values,k,h):
    peakScores = []
    peakScoresMax = []  #
    peakPositionsMax = []
    peakScoresMin = []
    peakPositionsMin = []
    #h = float(std_text.get("1.0", tk.END).strip('\n'))
    #k = int(window_size_average_text.get("1.0", tk.END).strip('\n'))

    for point in range(k,len(values)-k):
        value = values[point]
        signed_top = 0
        signed_bot = 0
        for i in range(k):
            signed_bot += value - values[point - i]
            signed_top += value - values[point + i]
        ps = (signed_top / k + signed_bot / k) / 2
        if (ps > 0):
            peakScoresMax.append(ps)
            peakPositionsMax.append(point)
        else:
            peakScoresMin.append(ps)
            peakPositionsMin.append(point)

    std_max = np.std(peakScoresMax, dtype=np.float64)
    mean_max = np.average(peakScoresMax)
    for i in range(len(peakScoresMax)):
        if (peakScoresMax[i] - mean_max > h * std_max):
            peakScores.append([peakPositionsMax[i], peakScoresMax[i]])

    std_min = np.std(peakScoresMin)
    mean_min = np.average(peakScoresMin)
    for i in range(len(peakScoresMin)):
        if (peakScoresMin[i] - mean_min < h * std_min * -1):
            peakScores.append([peakPositionsMin[i], peakScoresMin[i]])

    peakScores.sort(key=itemgetter(0))

    return peakScores

#Calculate a Peak Score for each pre-filtered local extrema. Peak Score is calculated by maximum distance around extrema
def calculatePeakScorePeaksMD(values,k,h,w):
    peakScores = []
    peakScoresMax = []  #
    peakPositionsMax = []
    peakScoresMin = []
    peakPositionsMin = []
    #filter_window_size = int(window_size_argrel_text.get("1.0", tk.END).strip('\n'))
    max_type = max_dropdown_variable.get()
    maximums = scipy.signal.argrelextrema(values, np.greater, 0, w)
    minimums = scipy.signal.argrelextrema(values, np.less, 0, w)
    #h = float(std_text.get("1.0", tk.END).strip('\n'))
    #k = int(window_size_average_text.get("1.0", tk.END).strip('\n'))

    if max_type == "All" or max_type == "Maximums":
        for point in maximums[0]:
            if point > k:
                if point < len(values) - k:
                    value = values[point]
                    signed_top = []
                    signed_bot = []
                    for i in range(k):
                        signed_bot.append(value - values[point - i])
                        signed_top.append(value - values[point + i])
                    ps = (max(signed_top) + max(signed_bot)) / 2
                    if (ps > 0):
                        peakScoresMax.append(ps)
                        peakPositionsMax.append(point)

        std_max = np.std(peakScoresMax, dtype=np.float64)
        mean_max = np.average(peakScoresMax)
        for i in range(len(peakScoresMax)):
            if (peakScoresMax[i] - mean_max > h * std_max):
                peakScores.append([peakPositionsMax[i], peakScoresMax[i]])

    if max_type == "All" or max_type == "Minimums":
        for point in minimums[0]:
            if point > k:
                if point < len(values) - k:
                    value = values[point]
                    signed_top = []
                    signed_bot = []
                    for i in range(k):
                        signed_bot.append(value - values[point - i])
                        signed_top.append(value - values[point + i])
                    ps = (min(signed_top) + min(signed_bot)) / 2
                    if (ps < 0):
                        peakScoresMin.append(ps)
                        peakPositionsMin.append(point)

        std_min = np.std(peakScoresMin)
        mean_min = np.average(peakScoresMin)

        for i in range(len(peakScoresMin)):
            if (peakScoresMin[i] - mean_min < h * std_min * -1):
                peakScores.append([peakPositionsMin[i], peakScoresMin[i]])

    peakScores.sort(key=itemgetter(0))

    return peakScores


#region tkinter global objects

root = Tk()
root.geometry('600x500')
file_label = tk.Label(root, text="")
loading_label = tk.Label(root, text="")

electrode_picker_label = tk.Label(root, text="Electrode Picker")
spike_detection_label = tk.Label(root, text="Spike Detection")
visualization_label = tk.Label(root, text="Visualization")
signature_matching_label = tk.Label(root, text="Signature Matching")
spacer1 = tk.Label(root, text="")
spacer2 = tk.Label(root, text="")
spacer3 = tk.Label(root, text="")

#region dropdown menus
x_dropdown_label = tk.Label(root, text="X")
x_dropdown_variable = tk.StringVar(root)
x_dropdown_variable.trace("w", refreshPeaks)
x_dropdown_choices = []
x_dropdown = tk.OptionMenu(root, x_dropdown_variable, x_dropdown_choices)
y_dropdown_label = tk.Label(root, text="Y")
y_dropdown_variable = tk.StringVar(root)
y_dropdown_variable.trace("w", refreshPeaks)
y_dropdown_choices = []
y_dropdown = tk.OptionMenu(root, y_dropdown_variable, y_dropdown_choices)
max_dropdown_label = tk.Label(root, text="Maxima")
max_dropdown_variable = tk.StringVar(root)
max_dropdown_choices = ['All','Maximums','Minimums']
max_dropdown_variable.set('All')
max_dropdown = tk.OptionMenu(root, max_dropdown_variable, *max_dropdown_choices)
method_dropdown_label = tk.Label(root, text="Spike Detection Method")
method_dropdown_variable = tk.StringVar(root)
method_dropdown_choices = ['Average Distance','Maximum Distance']
method_dropdown_variable.set('Average Distance')
method_dropdown = tk.OptionMenu(root, method_dropdown_variable, *method_dropdown_choices)

#region Text Inputs
window_size_average_label = tk.Label(root, text="Window size PeakScore")
window_size_average_text = tk.Text(root, height=1, width=10)
window_size_argrel_label = tk.Label(root, text="Window size Filter")
window_size_argrel_text = tk.Text(root, height=1, width=10)
std_label = tk.Label(root, text="std")
std_text = tk.Text(root, height=1, width=10)
lam_label = tk.Label(root, text="lambda")
lam_text = tk.Text(root, height=1, width=10)
threshold_big_label = tk.Label(root, text="Threshold PeakScore")
threshold_big_text = tk.Text(root, height=1, width=10)

#endregion

#region Buttons


button_calculate_peaks = tk.Button(root, text="Detect local Spikes", command=onCalculatePeaks)
button_calculate_all_peaks = tk.Button(root, text="Detect all Spikes", command=onCalculateAllPeaks)
button_show_electrode_plot = tk.Button(root, text="Show Electrode Plot", command=onShowElectrodePlot)
button_show_frequency_heatmap = tk.Button(root, text="Show Heatmap", command=onShowPeakFrequencyHeatmap)
button_create_excel_file = tk.Button(root, text="Create Excel File", command=onCreateExcelFile)
#button_perform_smd = tk.Button(root, text="Calculate all SMD", command=onCalculateSMD)
#button_perform_local_smd = tk.Button(root, text="Calculate local SMD", command=calculate_smd_one_electrode)

#endregion

#endregion

#region main

def startGUI():
    root.title("Spike Signature Matching")
    menubar = tk.Menu(root)

    filemenu = tk.Menu(menubar)
    filemenu.add_command(label="New cmcr", command=onFilemenuOpen)
    filemenu.add_command(label="Save", command=onSaveFile)
    filemenu.add_command(label="Load", command=onLoadFile)
    filemenu.add_command(label="Exit", command=root.quit)

    menubar.add_cascade(label="File", menu=filemenu)

    root.config(menu=menubar)

    rowcount = 0
    loading_label.grid(row=rowcount, column=0, columnspan=3)
    rowcount+=1

    file_label.grid(row=rowcount,column=0,columnspan=3)
    rowcount += 1

    electrode_picker_label.grid(row=rowcount,column=0,columnspan=3,padx=10)
    rowcount += 1

    x_dropdown.grid(row=rowcount, column=1)
    x_dropdown_label.grid(row=rowcount, column=0)
    rowcount += 1

    y_dropdown.grid(row=rowcount, column=1)
    y_dropdown_label.grid(row=rowcount, column=0)
    button_show_electrode_plot.grid(row=rowcount, column=2)
    rowcount += 1

    spacer1.grid(row=rowcount,column=0)
    rowcount += 1

    spike_detection_label.grid(row=rowcount,column=0,columnspan=3,padx=10)
    rowcount += 1

    method_dropdown.grid(row=rowcount,column=1)
    method_dropdown_label.grid(row=rowcount,column=0)
    rowcount += 1

    max_dropdown_label.grid(row=rowcount, column=0)
    max_dropdown.grid(row=rowcount, column=1)
    rowcount += 1

    window_size_argrel_text.grid(row=rowcount, column=1)
    window_size_argrel_text.insert(tk.INSERT, '25')
    window_size_argrel_label.grid(row=rowcount, column=0)
    rowcount += 1

    window_size_average_text.grid(row=rowcount, column=1)
    window_size_average_text.insert(tk.INSERT, '100')
    window_size_average_label.grid(row=rowcount, column=0)
    rowcount += 1

    std_text.grid(row=rowcount, column=1)
    std_text.insert(tk.INSERT, '2')
    std_label.grid(row=rowcount, column=0)
    button_calculate_peaks.grid(row=rowcount, column=2)
    button_calculate_all_peaks.grid(row=rowcount, column=3)
    rowcount += 1

    spacer2.grid(row=rowcount,column=0)
    rowcount += 1

    visualization_label.grid(row=rowcount,column=0,columnspan=3,padx=10)
    rowcount +=1

    threshold_big_text.grid(row=rowcount, column=1)
    threshold_big_text.insert(tk.INSERT, '70')
    threshold_big_label.grid(row=rowcount, column=0)
    button_show_frequency_heatmap.grid(row=rowcount, column=2)
    button_create_excel_file.grid(row=rowcount, column=3)
    rowcount += 1

    spacer3.grid(row=rowcount, column=0)
    rowcount += 1

    signature_matching_label.grid(row=rowcount,column=0,columnspan=3,padx=10)
    rowcount += 1

    lam_text.grid(row=rowcount, column=1)
    lam_text.insert(tk.INSERT, '0.9')
    lam_label.grid(row=rowcount, column=0)
    #button_perform_smd.grid(row=rowcount, column=3)
    #button_perform_local_smd.grid(row=rowcount, column=2)
    rowcount += 1



    root.mainloop()

def showElectrodePlot(electrode):
    global td
    tr = td.tickrate
    xs = range(len(electrode.values))
    xs = xs*tr
    xs = xs / 1000000
    fig, ax = plt.subplots()
    markers_on = electrode.getPeakPositions()
    ax.plot(xs, electrode.values,'-o', markevery=markers_on, color='tab:blue', mfc='r')
    ax.set_xlabel('time in s')
    ax.set_ylabel('μV')
    ax.set_title('X: ' + str(electrode.x) + ' Y: ' + str(electrode.y))
    plt.show()

if __name__ == '__main__':


    startGUI()
