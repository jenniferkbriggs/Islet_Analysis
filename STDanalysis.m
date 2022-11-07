function CellTC = STDanalysis(images, CellMask,opts)
%Jennifer Briggs 02.2022
%This function takes a time course of pixels with masks differentiating
%cells and removes pixels within that mask that do not correlate with the
%rest of the pixels (e.g. likely not a part of the cell.)

%images is the image file in matrix form (x,y,t). Currently only accepting
%3D matrices (e.g. no z stacks) so z stacks should be fed in separately. 

%opts is a structure of options:
%opts.figs = 1 if want a gif to see the pixels be removed

    CellMasksave = CellMask; %save old cell mask
    
    numcells = max(max(CellMask))
    for i = 1:numcells
        cctt = 1;
        badpix = [2 2]; %arbitrary array that will be filled with the bad pixel locations
        disp(i)
        while ~isempty(badpix)
        TCMask = CellMask; %Pulls in CellMask array
        TCMask(TCMask ~= i) = 0; %Gets rid of all masks besides current one
        MaskedIMGstack = images.*logical(TCMask); %Applies current mask to all frames of timecourse
        
        %get rid of for loop
        [x,y] = find(TCMask);
         %Middleimage = MaskedIMGstack(round(mean(x))-5:round(mean(x))+5, round(mean(y))-5:round(mean(y))+5, :);
        [rr, cc] = find(mean(MaskedIMGstack,3));
        for ii = 1:size(MaskedIMGstack,3)
            TCnoZero = MaskedIMGstack(:,:, ii); %Pulls in current frame
            TCnoZero = TCnoZero(TCnoZero>0); %Accounts for any zeros from preallocation
            TCcheck(ii,:) = TCnoZero;
            %TCmiddle(ii) = mean2(Middleimage(:,:,ii));
            TC(ii) = mean(TCnoZero);
        end
        
        TCchecksm = smoothdata(TCcheck); %Smooth data
        if size(TCchecksm,2) >0
        for c=1:size(TCchecksm,2)
                [C1] = xcorr(TCchecksm(:,c),median(TCchecksm'),3);
                Maxcor(c) =  max(C1);       
        end
        %pixel is bad if it shows no oscillations
        
        %badpix = find(Maxcor<5e6);
        %if length(badpix) > .75*size(TCcheck,2)       
               badpix = find(Maxcor<mean(Maxcor) - 1.5*std(Maxcor))
               cctt = cctt +1;
       %end
       if cctt < 8  

        figure, title(['Itteration' num2str(cctt)]),
        plot(rr,cc,'o'), 
        %saveas(gcf, ['Itteration' num2str(cctt), '.png']) 
        for bb= 1:length(badpix)
        hold on, plot(rr(badpix(bb)),cc(badpix(bb)),'ro','markerfacecolor','r')
        CellMask(rr(badpix(bb)),cc(badpix(bb))) = 0;
        end
       else
    %        disp('Bad Pixels are too many')
    %        disp(savepath)
    %        fig = figure, plot(TCchecksm)
    %        x = input('Keep deleting or move on? 1 for move on, 0 for keep deleting')
    %        close(fig)
    %         if x == 0
    %             cctt = 1;
    %         else
    %             cctt = 5;
                badpix = [];
    %         end
       end
            clear TCcheck Maxcor Maxcor2

            %figure, plot(TCchecksm)
         else
           disp('No pixels')
           badpix = [];
           clear TCcheck
        end
        end
close all 

        CellTC(:,i) = TC; %Updates CellTC array with mean intensities from each frame of each cell
        PlotLabels{i} = ['Cell' num2str(i)]; %Updates labels with current cell number for legend
    end
end