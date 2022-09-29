% Run Network analysis: 

%load csv: 
capath = '/Users/levittcl/Documents/Hub Analysis/hubanalysistest.mat\CaWaveForm.mat'%%%% Add .mat path here
ca = importdata(capath);

    a = find(kperc(l).Islet(i).data  > 60);
    hubthreshold = (a(1));