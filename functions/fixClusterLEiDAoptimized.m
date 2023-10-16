function [PTRANSITION,Pstates,transferMatrix] = fixClusterLEiDAoptimized(BOLD,numClusters,Center,TR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Jakub Vohryzek 06/11/2018
% from Gustavo Deco awakening code
% for the current project: Treatment-Resistant Depression in psilocybin treatment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[numAreas, tMax] = size(BOLD);

% Preallocate variables to save FC patterns and associated information
leadingEigV1 = zeros(tMax,numAreas); % First leading eigenvectors
% % % leadingEigV1V2 = zeros(tMax,2*numAreas); % All leading eigenvectors

%% BANDPASS FILTER SETTINGS

fnq = 1/(2*TR);                 % Nyquist frequency
flp = 0.04; %0.02; % 0.04;                    % lowpass frequency of filter (Hz)
fhi = 0.07; %0.1; % 0.07;                    % highpass
Wn = [flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k = 2;                          % 2nd order butterworth filter
[bfilt,afilt] = butter(k,Wn);   % construct the filter
clear fnq flp fhi Wn k

%% CODE
tAll = 0; % Index of time (starts at 0 and will be updated until n_Sub*Tmax)
phaseBOLD = zeros(numAreas,tMax);

%% Get the BOLD phase using the Hilbert transform
for seed=1:numAreas
    BOLD(seed,:)= BOLD(seed,:)-mean(BOLD(seed,:));
    signalFilt =filtfilt(bfilt,afilt,BOLD(seed,:));
    phaseBOLD(seed,:) = angle(hilbert(signalFilt));
end

%% Calculate the Instantaneous FC (BOLD Phase Synchrony)
for t=1:tMax
    iFC=zeros(numAreas);
    for n=1:numAreas
        for p=1:numAreas
            iFC(n,p) = cos(phaseBOLD(n,t)-phaseBOLD(p,t));
        end
    end
    
    % Get the leading eigenvector
    % [V1 ~]=eigs(iFC,1);

    [VV DD]=eigs(iFC,2);
%     CONSISTENCY with LEiDA and thus not weighted by the eigenvalue
%     d1=DD(1,1)/sum(diag(DD));
%     d2=DD(2,2)/sum(diag(DD));
%     V1=d1*VV(:,1);
%     V2=d2*VV(:,2);
    V1 = VV(:,1);
% % %     V2 = VV(:,2);

    % Make sure the largest component is negative
    if mean(V1>0)>.5
        V1=-V1;
    elseif mean(V1>0)==.5 && sum(V1(V1>0))>-sum(V1(V1<0))
        V1=-V1;
    end
    
% % %     if mean(V2>0)>.5
% % %         V2=-V2;
% % %     elseif mean(V2>0)==.5 && sum(V2(V2>0))>-sum(V2(V2<0))
% % %         V2=-V2;
% % %     end
    
    % Save V1 from all frames in all fMRI sessions in Leading eig
    tAll=tAll+1; % Update time
    leadingEigV1(tAll,:) = V1;
    % leadingEigV1V2(tAll,:) = vertcat(V1,V2);
end

clear iFC V1 Phase_BOLD

%% CLUSTERING THE LEADING EIGENVECTOR(S)
IDX = zeros(tAll,1);

for t=1:tAll
    for j = 1:numClusters
        % JAKUB
     di(j) = sqrt(sum((leadingEigV1(t,:)-Center(j,:)).^2)); % comparing clusters to the empirical data
     % the ^2 is done after the sum - not exactly a euclidean distance -
     % does it make sense?
     % DECO
% % %      diDeco(j)=sqrt(sum(leadingEigV1(t,:)-Center(j,:)).^2); % comparing clusters to the empirical data
    end
    % JAKUB
    [aux, indmin] = min(di);
    IDX(t) = indmin;
    % DECO
% % %     [auxDeco, indminDeco] = min(diDeco);
% % %     IDXDeco(t) = indminDeco;
end

% JAKUB OPTION 1: Fractional Occupancy
Pstates=zeros(1,numClusters);
for c = 1:numClusters
    Pstates(c) = mean(IDX == c);
end
Pstates=Pstates/sum(Pstates);

% DECO OPTION 1: Fractional Occupancy
% % % PstatesDeco=zeros(1,numClusters);
% % % for c = 1:numClusters
% % %     PstatesDeco(c) = mean(IDXDeco == c);
% % % end
% % % PstatesDeco = PstatesDeco/sum(PstatesDeco);

% OPTION 1
PTRANSITION=zeros(numClusters,numClusters);
i=1;
for c1=1:numClusters
    j=1;
    for c2=1:numClusters
        sumatr=0;
        for t=1:length(IDX)-1
            if IDX(t)==c1 && IDX(t+1)==c2
                sumatr=sumatr+1;
            end
        end
        if length(find(IDX(1:length(IDX)-1)==c1)) ~= 0
            PTRANSITION(i,j)=sumatr/(length(find(IDX(1:length(IDX)-1)==c1))); % normalised by the rows
        end
        j=j+1;
    end
    i=i+1;
end
% OPTION 2
transferMatrix = zeros(numClusters,numClusters);
for tp = 2:length(IDX)
    transferMatrix(IDX(tp-1),IDX(tp)) = transferMatrix(IDX(tp-1),IDX(tp)) + 1;
end
transferMatrixAllnorm = transferMatrix/sum(sum(transferMatrix)); % normalised by the whole matrix
transferMatrix = transferMatrix./sum(transferMatrix,2); % normalised by the rows as in Gus's code