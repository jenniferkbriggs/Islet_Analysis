% Run Network analysis: 

%load cvs: 
capath = %%%% Add .csv path here
ca = importdata(capath);

    a = find(kperc(l).Islet(i).data  > 60);
    hubthreshold = (a(1));