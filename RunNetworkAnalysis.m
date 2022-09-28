% Run Network analysis: 
% Jennifer Briggs 2022


%load cvs: 
capath = %%%% Add .csv path here
ca = importdata(capath);

Threshold = 0.99 %Here you put the correlation threshold to draw an edge
Opts.figs = 1 %Set 1 if you want figures, 0 if not


[N, Adj, kpercent, histArrayPercShort,pval,Rij,s] = NetworkAnalysis(ca, Threshold, Opts)

a = find(kperc > 60); hubthreshold = (a(1)); %Find degree threshold for hubs
degree = sum(Adj); %find degree of each cell
Hubs = find(degree>hubthreshold) %List each hug