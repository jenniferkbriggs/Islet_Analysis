%% Show network dot stick
close all 
clear all
clc
addpath('~/GitHub/UniversalCode');

% THINGS YOU CHANGE
filepath = %input directory where CaWaveForm.mat is 
imagepath = %input directory where Imaging.mat is
savename = %input where to save the data
Thr = .9   %input correaltion threshold for network analysis
TitleChoice = 'input title choice here' %Input title choice here
%Input start and end time of video
starttime = Vidinfo(illy).starttime(lg); 
endtime = Vidinfo(illy).endtime(lg);

% Load things: 
load([filepath '\CaWaveForm.mat'])
load([imagepath '\' 'Imaging.mat'])
load([filepath '\' 'Masks.mat'])       %load masks
load([imagepath '\' 'CellNumber.mat'])

% 
zstacks = 3         %how many z stacks
zz = 2              %where to start on z stack
cachannel = 3       %where is the calcium channel
howmanychannel = 3  %how many imaging channels


% Load video
pics=R{1};
pics=pics(:,1);

for i=1:length(pics)
    IMG(:,:,i)=pics{i};
end
pn = length(pics);

for i=1:pn
    IMG(:,:,i)=pics{i};
end

try
    for i=1:length(pics)
        T(i)=R{4}.getPlaneDeltaT(0, i-1).value;
    end
catch
    T=0:0.5:pn*0.5;
end
T = double(T);
T = T(cachannel:howmanychannel:end);
T = T(1:zstacks:end);

if starttime == -1
    st=1;
else
    st = starttime;
end

if endtime == -1
    ed=length(T);
else
    ed=endtime;
end

T = T(st:ed);
images=double(IMG); % converts images to double precision
images = images(:,:,cachannel:howmanychannel:end);
RawImg=images(:,:,1); % assigns the first frame of the video to RawImg variable
images = images(:,:,zz:zstacks:end);
images = images(:,:,st:ed-1);

sx=size(images,1);
sy=size(images,2);
sz=length(T);
for i=1:size(images,3)
    images(:,:,i)=medfilt2(images(:,:,i),[5 5]); %applies filter to clean up images
end
toc
end
ImAv = sum(images,3); %compresses all frames into single array of intensities
HSV = ones(sx,sy,3); %preallocates a 3 dimensional array
ImAvn = ImAv/max(ImAv(:));
HSV(:,:,3) = ImAvn.^0.8; %evens out intensity across images
HSV(:,:,1) = 0.3333;%converts image to green image
RGB2 = hsv2rgb(HSV); %converts to rgb image

%%
for ll = 1:numcells
  
    
    [xp,yp] = find(CellMask == ll);
    x(ll) = mean(xp);
    y(ll) = mean(yp);
end

bad = find(isnan(x))

x(bad) = [];
y(bad) = [];
CellTC(:, bad) = [];
numcells = numcells - length(bad)
%%
fig = figure
imshow(RGB2)%, hold on

Adj = corr(CellTC);
Adj = Adj > Thr;
Adj = Adj - diag(diag(Adj));
AdjacencyGraph = graph(Adj);
Conn = sum(Adj);

Conarray = [0:max(Conn)];
Conarray2 = Conarray/length(Conarray)*100;
try
    hubth = Conarray(Conarray2 > 60); hubth = min(hubth); %hub threshold
    Hubby = find(Conn >= hubth); %Find Hubs: Number of cells > 60% of Islet Links
    Nodec = repmat([.175 .54 .60],numcells,1); %color of nodes
    for lll = 1:length(Hubby) %load nodes 
        Nodec(Hubby(lll),:)=[.98 .122 .157];
    end
catch
    Nodec = repmat([.175 .54 .60],numcells,1);
end

fig2 = figure
p = plot(AdjacencyGraph, 'Xdata',y,'YData',x, 'EdgeColor', 'b', 'NodeColor',Nodec,'MarkerSize',8, 'LineWidth',1 )
p.NodeLabel = [];
set(gca, 'YDir','reverse')
title([TitleChoice])



set(gca, 'Visible', 'off')
set(gca, 'xticklabels', [])
set(gca, 'yticklabels', [])
% 
set(fig, 'Position', [100 100 1000 800])
set(fig2, 'Position', [100 100 1000 800])


saveas(fig, [savename '.png'])
saveas(fig2, [savename '.png'])

clearvars x y
close(fig)
close(fig2)


