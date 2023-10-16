function [h,p,ktest] = FCDhistktest(groupA,groupB,numBin)
% This functions calculates the h, pvalue and kvalue between two
% distribution with varying binSize
    hPRE = histogram(groupA,numBin,'Normalization','Probability');
    ValuePRE = hPRE.Values;
    hPOST = histogram(groupB,numBin,'Normalization','Probability');
    ValuePOST = hPOST.Values;
    [h,p,ktest] = kstest2(ValuePRE,ValuePOST);

end

