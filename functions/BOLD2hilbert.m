function [ BOLD_processed, Phase_BOLD, staticFC] = BOLD2hilbert( BOLD, numAreas, bfilt, afilt, excTp )

%   Detailed explanation goes here
        for seed = 1:numAreas
            tmp = demean(BOLD(seed,:)); % demean(detrend(signalData(seed,:))); maybe as Deco's code?
            signalFilt = filtfilt(bfilt,afilt,tmp);
            BOLD_processed(seed,:) = signalFilt((excTp+1):(end-excTp)); % excluding first and last 3 tps

            % BOLD_processed(seed,:) = signalFilt; % excluding first 3 tps

            Phase_BOLD(seed,:) = angle(hilbert(signalFilt((excTp+1):(end-excTp)))); % excluding first and last 3 tps
        end     
        staticFC = corrcoef(BOLD_processed');

end

