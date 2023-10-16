function [ ts_processed, Phase_BOLD, staticFC] = ts2filtered( BOLD, numAreas, bfilt, afilt, excTp )

%   Detailed explanation goes here
        for seed = 1:numAreas
            tmp = detrend(BOLD(seed,:)-nanmean(BOLD(seed,:)));
            signalFilt = filtfilt(bfilt,afilt,tmp);
            ts_processed(seed,:) = signalFilt((excTp+1):(end-excTp)); % excluding first and last 3 tps

            % BOLD_processed(seed,:) = signalFilt; % excluding first 3 tps

            Phase_BOLD(seed,:) = angle(hilbert(signalFilt((excTp+1):(end-excTp)))); % excluding first and last 3 tps
        end     
        staticFC = corrcoef(ts_processed');

end

