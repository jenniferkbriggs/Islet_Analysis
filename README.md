# Islet_Analysis
This repository contains cleaned and commented code for general analyses conducted on islet timeseries in the Benninger Lab, University of Colorado Anschutz Medical Campus.

## Loading data
The function file **LoadData** is an easy way to load time series data and get it into a csv form ready for the rest of the analyses. Currently it can load .csv and .cmcr files. 

### Dependencies: 
easygui, McsPy

## Network analysis in Python
### Dependencies: 
The code depends on the following python packages: csv, collections, networkx, numpy, matplotlib.pyplot, scipy, easygui, pandas, math, tkinter, pyvis. **If you get an error with that looks like: Exception has occurred: ModuleNotFoundError: no module named '<modulename>'. Then that module must be downloaded using either Conda or pip. To download the module, go to the terminal and type (for example) 'conda install -c conda-forge easygui'.

### To run a general network analysis, you will use the code **Run_Network.py**. 
The first section (lines 5-17) is title options for you to change. This is wheter you define whether you would like figures, where to save the file to.




