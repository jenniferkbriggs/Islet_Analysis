%% THIS PROGRAM IMPORTS CALCIUM IMAGING FILES FOR CELL-BY-CELL ANALYSIS OF CALCIUM TRACE TO IDENTIFY Cell Network
%% REAL-TIME USER INPUT REQUIRED
%% Jennifer Briggs, Feb 2021
%% HOUSEKEEPING
close all
clear all
clc


cachannel = 1;
nuchannel = 2;
howmanychannel = 1;

filename = '/Users/levittcl/Documents/Research/DATA/gCAMP - Calcium Imaging/2022_06_16/gCAMP_islet1_control_2mM_11mM_KCl.czi' %file name can direct to '.mat' analysis file or imaging file
savepath = '/Users/levittcl/Documents/Hub Analysis/'
savename = '/hubanalysistest.mat'



numad = 1;
starttime = -1
endtime = -1
gstart = 1

    
%% LOADING THE CA IMAGE FILE
if exist('R') == 0
try
    load([filename], 'R', 'zz', 'zstacks', 'cachannel', 'howmanychannel')
catch
    R = bfopen([filename]); % Uses bfopen program to open .czi/.lsm image files
    try
        save([savepath savename],  '-V7.3')
    catch
        mkdir([savepath])
        save([savepath savename],  '-V7.3')
    end
 end
end


%% 
% omeMeta = 0
{1,4};
% voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); % in µm
% voxelSizeXdouble = voxelSizeX.doubleValue();                                  % The numeric value represented by this object after conversion to type double                             
tic
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


if starttime == -1
    st = 1;
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

clear pics R IMG;
% clear omeMeta;

output_dir = savepath;
toc

%% DECLARING IMAGE PROPERTIES
tic
images = images(:,:,st:ed-1);
sx=size(images,1);
sy=size(images,2);
sz=length(T)-1;
for i=1:sz
    images(:,:,i)=medfilt2(images(:,:,i),[5 5]); %applies filter to clean up images
end

%Makes all frames 1fps
% fps =Vidinfo(g).fps(ff);
fps = 1;
 ct =1
for i = 1:fps:sz
    images(:,:,ct) = images(:,:,i);
    ct = ct+1;
end
images = images(:,:,1:ct-1);
%ImAv = sum(images,3); %compresses all frames into single array of intensities
ImAv = mean(images,3); 
HSV = ones(sx,sy,3); %preallocates a 3 dimensional array
ImAvn = ImAv/max(ImAv(:));
HSV(:,:,3) = ImAvn.^0.8; %evens out intensity across images
HSV(:,:,1) = 0.3333;%converts image to green image

RGB2 = hsv2rgb(HSV); %converts to rgb image

RGB= im2gray(RGB2);
RGB = imadjust(RGB);
OGFig = figure(1);
imshow(RGB)
toc
%% 

%%
%% MASKING ISLET FOR CELLS



%% User Draws ROIs around each cell within the islet
% ROIs are saved in "CellMask" array and called back to throughout analysis
CellMask = double(zeros(sx,sy));
numcells = 1;

% 
try %try to load in cell masks if already in loadpath
    try
    load([loadpath '\Masksv2.mat'])
    catch
    load([loadpath '\Masks.mat'])
    end
    load([loadpath '\CellNumber.mat'])
    
    ImAv = sum(images.*logical(CellMask),3); %compresses all frames into single array of intensities
    HSV = ones(sx,sy,3); %preallocates a 3 dimensional array
    ImAvn = ImAv/max(ImAv(:));
    HSV(:,:,3) = ImAvn.^0.8; %evens out intensity across images
    HSV(:,:,1) = 0.3333;%converts image to green image
    RGB2 = hsv2rgb(HSV); %converts to rgb image
    imshow(RGB2)   

catch
    NoSigFig = figure('Name','Draw ROIs Around Cells of Interest');
    imshow(RGB);
    keyboard
    k = 1;
    while k > 0
        disp('Draw ROIs Around Cells of Interest')
        ROIMask = imfreehand(); %User draws region around cell
        ROIMask = createMask(ROIMask); %Mask is created from drawn region
        CellMask = CellMask + ROIMask.*numcells; %CellMask array is updated with new mask; new mask is multiplied by the cell label before updating
        CellMask(find(CellMask>numcells)) = numcells; %If a region is overlapped, it is instead attributed to the most recent region
        UInp = input('Select Additonal ROI? 1=Yes, 0=No, 2=Redraw \n'); % No preallocated number of ROIs, User input required to continue ROI drawing or to stop
        if UInp == 1 %User input to determine if another region is drawn
            k = k+1;
            numcells = numcells+1;
        elseif UInp == 2
            numcells = numcells;
            k = k;
        elseif UInp ~= 1
            k = 0;
        end
    end
    save([savepath '\Masks.mat'],'CellMask')
    save([savepath '\CellNumber.mat'],'numcells')
    close(NoSigFig);
    clear NoSigFig;  


    %run STD analysis and extract calcium wave form
    CellTC = STDanalysis(images, CellMask);

    
    % Plotting traces for entire time course
    TCFig = figure('Name','Average Intensity Over Time');
    plot(CellTC);
    
    %clear images MaskedIMGstack;
    try
    saveas(TCFig,[savepath '/Cellintestiy.tif']); %Saves figure of each cell's timecourse
    save([savepath savename '/CaWaveForm.mat'],'CellTC')
    catch
        mkdir(savepath)
        saveas(TCFig,[savepath savename '/Cellintestiy.tif']); %Saves figure of each cell's timecourse
    end
    save([savepath savename '/CaWaveForm.mat'],'CellTC')
   
end


