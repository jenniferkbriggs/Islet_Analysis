% Run Network analysis: 
% Jennifer Briggs 2022
close all
clear all

%load cvs: 
capath = '/Users/levittcl/Documents/Research/Projects/SST/Calcium Delta Cells/3_12_2023/islet 4 5mM 11mM/CaWaveForm.mat' %%%% Add .mat path here
ca = importdata(capath);

%To set the threshold, either manually set: 

%1) ----------
Threshold = 0.9 %Here you put the correlation threshold to draw an edge

%Or determine the threshold using the average degree or a threshold that
%gives you a scale free degree distribution

%2) --------------
Opts.Method = 'Degree'
Opts.avDeg = 6; %Set the average degree here

Opts.Method = 'Scale-Free' 
%%Set the bounds for max and min average degree for scale free: 
Opts.Max = 20
Opts.Min = 2
%NOW RUN
Threshold = findoptRth(calcium, Opts)


Opts.figs = 0 %Set 1 if you want figures, 0 if not

[degree, Adj, kpercent, histArrayPercShort,pval,Rij,s] = NetworkAnalysis(ca, Threshold, Opts)

a = find(kpercent > 60); hubthreshold = (a(1)); %Find degree threshold for hubs
Hubs = find(degree>hubthreshold) %List each hug <3