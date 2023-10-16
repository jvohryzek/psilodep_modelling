function [PTRANSITION,Pstates] = LEiDA_fix_cluster(BOLD,NumClusters,Center,TR)

[numAreas, tMax] = size(BOLD);

% Preallocate variables to save FC patterns and associated information
leadingEig = zeros(tMax,2*numAreas); % All leading eigenvectors

% Bandpass filter settings
fnq = 1/(2*TR);                 % Nyquist frequency
flp = 0.02 % 0.04;                    % lowpass frequency of filter (Hz)
fhi = 0.1 % 0.07;                    % highpass
Wn = [flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k = 2;                          % 2nd order butterworth filter
[bfilt,afilt] = butter(k,Wn);   % construct the filter
clear fnq flp fhi Wn k

t_all = 0; % Index of time (starts at 0 and will be updated until n_Sub*Tmax)

phaseBOLD = zeros(numAreas,tMax);

% Get the BOLD phase using the Hilbert transform
for seed=1:numAreas
    BOLD(seed,:)=BOLD(seed,:)-mean(BOLD(seed,:));
    signalFilt =filtfilt(bfilt,afilt,BOLD(seed,:));
    phaseBOLD(seed,:) = angle(hilbert(signalFilt));
end

for t=1:tMax
    
    %Calculate the Instantaneous FC (BOLD Phase Synchrony)
    iFC=zeros(numAreas);
    for n=1:numAreas
        for p=1:numAreas
            iFC(n,p) = cos(phaseBOLD(n,t)-phaseBOLD(p,t));
        end
    end
    
    % Get the leading eigenvector
%      [V1 ~]=eigs(iFC,1);

    [VV, DD]=eigs(iFC,2);
    d1=DD(1,1)/sum(diag(DD));
    d2=DD(2,2)/sum(diag(DD));
    V1=d1*VV(:,1);
    V2=d2*VV(:,2);
    % Make sure the largest component is negative
    if mean(V1>0)>.5
        V1=-V1;
    elseif mean(V1>0)==.5 && sum(V1(V1>0))>-sum(V1(V1<0))
        V1=-V1;
    end

    if mean(V2>0)>.5
        V2=-V2;
    elseif mean(V2>0)==.5 && sum(V2(V2>0))>-sum(V2(V2<0))
        V2=-V2;
    end
    
    % Save V1 from all frames in all fMRI sessions in Leading eig
    t_all=t_all+1; % Update time
    leadingEig(t_all,:)=vertcat(V1,V2);
end

clear signal_filt iFC V1 Phase_BOLD

%% 2 - Cluster the Leading Eigenvectors
IDX = zeros(t_all,1);

for t=1:t_all
    for j=1:NumClusters
        di(j)=sqrt(sum(leadingEig(t,:)-Center(j,:)).^2); % comparing clusters to the empirical data
    end
    [aux indmin]=min(di);
    IDX(t)=indmin;
end

Pstates=zeros(1,NumClusters);
for c=1:NumClusters
    Pstates(c)=mean(IDX==c);
end

Pstates=Pstates/sum(Pstates);

PTRANSITION=zeros(NumClusters,NumClusters);
i=1;
for c1=1:NumClusters
    j=1;
    for c2=1:NumClusters
        sumatr=0;
        for t=1:length(IDX)-1
            if IDX(t)==c1 && IDX(t+1)==c2
                sumatr=sumatr+1;
            end
        end
        if length(find(IDX(1:length(IDX)-1)==c1)) ~= 0
            PTRANSITION(i,j)=sumatr/(length(find(IDX(1:length(IDX)-1)==c1)));
        end
        j=j+1;
    end
    i=i+1;
end

 