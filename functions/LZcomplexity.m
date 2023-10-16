function [C_norm ,C] = LZcomplexity(BOLD_processed, numAreas)
        
% Calculate Lempel-Ziv Complexity
for i=1:numAreas
    tmp = BOLD_processed(i,:);
    C = LZ76(tmp>0); % setting all values above 0 to 1, else to 0 see Mediano et al. 2020 for discussion
    C_norm(i) = C*log2(size(tmp,2))/size(tmp,2);
end


end

