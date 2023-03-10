# Islet_Analysis
This repository contains cleaned and commented code for general analyses conducted on islet timeseries in the Benninger Lab, University of Colorado Anschutz Medical Campus.

Note: The python code is executed under virtual environment for version control purposes. The best way to avoid bugs is to run it in the same virtual environment.
The first time you clone this code, run: *conda env create -f environment.yml* in the terminal.
Then, before every run type *conda activate ./envs* to activate the environment before you run the files.

## Image Processing in MatLab
*This description is not complete*
The matlab file: *Extracting_individual_cells.m* is a script used to load calcium images and manually circle the cells. After all of the cells have been identified, the script will run *STD_analysis.m* which removes pixels from cells whose timecourses are too different from the rest of the cell. The script will output the final cell mask and the calcium csv file which contains timecourse for each cell.

## Network Analysis in Matlab
After you run *Extracting_individual_cells.m*, you will have a .mat file with all CaWaveForm.mat. To run general network analysis, you will open *RunNetworkAnalysis.m* and add the file location to this .mat file under capath

If you want to determine the optimal threshold for getting a log-log plot, input calcium into *findoptRth*. Please note, if you want to have the threshold based on number of links, please email me and I will add this to the code!

## Network analysis in Python
### General Network Analysis: **Run_Network.py**.
The first section (lines 5-17) is title options for you to change. This is wheter you define whether you would like figures, where to save the file to. ]

Try typing python3 Run_Network.py in the terminal and you should have success running this code.


### Dependencies:
The code depends on the following python packages: csvkit, collection, networkx, numpy, matplotlib.pyplot, scipy, easygui, pandas, python-math, tk, pyvis. **If you get an error with that looks like: Exception has occurred: ModuleNotFoundError: no module named '<modulename>'. Then that module must be downloaded using either Conda or pip. To download the module, go to the terminal and type (for example) 'conda install -c conda-forge easygui'.



## Signal Processing in Python
###  Filtering data with convolutional window: **Run_SignalProcessing.py**
This code will first load data using the dependency **LoadData.py** located in the git. The program will create a hamming window with size *windowSize*. *windowSize* is what you will change to alter the frequencies of the filter. Note that the window size is based on enumerated time, not real time. (e.g. the size 40 is 40 data points not 40 seconds.)

### Loading data
The function file **LoadData** is an easy way to load time series data and get it into a csv form ready for the rest of the analyses. Currently it can load .csv and .cmcr files.

### Dependencies:
easygui, McsPyDataTools
