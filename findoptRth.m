function [Thr] = findoptRth(calcium)
%This code allows you to find the optimal threshold for a network that both
%optimized fit to a scalefree graph and constraints for the average number of
%connections per islet
%Input: calcium
%Output: Threshold
% Jennifer Briggs, 2021

upper = 1 %set average number of connections

addpath('~/Documents/GitHub/UniversalCode/FMINSEARCHBND')
[Rij, pval]=corr(calcium);

fun = @(x)Lfunc(x,calcium,Rij); % Defines the cost function
numcells = size(calcium,2); 
Rij = Rij-diag(diag(Rij)); %removes ones on diagonal
Rijs = sort(Rij);
av = mean(Rijs(end-upper-1:end,:)); %calculates the average correlation coeffiecient for top 9 highest correlated cell pairs [9 x numcells]
avs = sort(av); %Gives average correlation for top 9 highest correlated cellpairs
ub = mean(avs(end-upper-1:end)); %Upper bound is set to be the average correlation for the top 9

x0 =  mean(avs(end-round(numcells/3):end)); %originial threshold


%run optimization: fminsearchbnd can be downloaded here https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon
Thr = fminsearchbnd(fun,x0, [0.6], ub); %Find threshold that most fits power law while staying within upper and lower bounds


end




function err = Lfunc(Threshold,calciumT,Rij)

numcells = size(calciumT,2);

Adj = Rij;
Adj(Adj >= Threshold) = 1;
Adj(Adj < Threshold) = 0;
% 
 Adj = Adj - diag(diag(Adj));        
 if mean(nonzeros(sum(Adj))) > 5 %here is where the lower bound is set!
%% 4. Determine number of "links" based on cov threshold
for i=1:numcells
    N (i,1) = nnz(Adj(:,i));  % N is matrix containing # of links for each cell (nnz returns number of nonzero elements)
end
% 


% %% 5. Creating a "probability of a cell to have k links" matrix

histArray=zeros(size(N))'; % prealocate
% a forloop to count how many times you encounter in a particular value:
for n=1:length(N)
    histArray(1,N(n)+1)=histArray(1,N(n)+1)+1; % every time you meet the particular value, you add 1 into to corresponding bin
end


histArrayPerc=histArray.*100/sum(histArray); % converting into % via dividing # of cell with specific # of links by total # of links 

m=find(histArray);    % returns all indexes of non-zero elements of the array
maxNonZeroLinks=max(m);   % number of links to consider when normalizing probabilty
k=1:1:maxNonZeroLinks;            % index of # of links (starting with 0 links)
kpercent=k.*100/(maxNonZeroLinks);  % convert # of links into % of limks
histArrayPercShort = histArrayPerc(1:maxNonZeroLinks);   % cropping the hisArray after the last non-0 value

histArrayShort = histArray(1:maxNonZeroLinks);
loghist = log(histArrayShort);
xlab = [1:length(histArrayPercShort)];
xlab(isinf(loghist))=[];
loghist(isinf(loghist))=[];


[s] = corrcoef(log(xlab),loghist)

err = abs(-1-s(2))
 else
     err = 2;
 end
end