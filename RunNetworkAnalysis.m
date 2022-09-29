% Run Network analysis: 
% Jennifer Briggs 2022
close all
clear all

%load cvs: 
capath = '/Users/levittcl/Documents/Hub Analysis/hubanalysistest.mat\CaWaveForm.mat' %%%% Add .mat path here
ca = importdata(capath);

Threshold = 0.9 %Here you put the correlation threshold to draw an edge
Opts.figs = 1 %Set 1 if you want figures, 0 if not


[degree, Adj, kpercent, histArrayPercShort,pval,Rij,s] = NetworkAnalysis(ca, Threshold, Opts)

a = find(kpercent > 60); hubthreshold = (a(1)); %Find degree threshold for hubs
Hubs = find(degree>hubthreshold) %List each hug <3