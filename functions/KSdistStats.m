function [statsKS] = KSdistStats(data)
% Calculates KSdist statistics as a function of bins

    % x*normal+y distributions

    FCDvaluesBefore = [data.remission.beforeRest.histFCD; data.no_remission.beforeRest.histFCD];
    FCDvaluesAfter = [data.remission.afterRest.histFCD; data.no_remission.afterRest.histFCD];
    statsKS = struct;
    ind =0;

    for i=10:100:10000
        ind = ind+1;
        % real
        [statsKS.hreal(ind),statsKS.preal(ind),statsKS.ktestreal(ind)] = FCDhistktest(FCDvaluesBefore,FCDvaluesAfter,i);
        % normal distribution
        [statsKS.hrand(ind),statsKS.prand(ind),statsKS.ktestrand(ind)] = FCDhistktest(randn(size(FCDvaluesBefore)),randn(size(FCDvaluesAfter)),i);
        % real with different mean and std. dev
        a=5;b=500;
        [statsKS.hrandshifted(ind),statsKS.prandshifted(ind),statsKS.ktestrandshifted(ind)] = FCDhistktest(a.*randn(size(FCDvaluesBefore'))+b,a.*randn(size(FCDvaluesAfter))+b,i);
        % real after responders and non-responders
        [statsKS.hrealrem(ind),statsKS.prealrem(ind),statsKS.ktestrealrem(ind)] = FCDhistktest(data.remission.afterRest.histFCD,data.no_remission.afterRest.histFCD,i);
        % real before responders and non-responders against after responders
        [statsKS.hrealrembfvsafr(ind),statsKS.prealrembfvsafr(ind),statsKS.ktestrealrembfvsafr(ind)] = FCDhistktest(FCDvaluesBefore,data.remission.afterRest.histFCD,i);
        % real before responders and non-responders against after non-responders
        [statsKS.hrealrembfvsafnr(ind),statsKS.prealrembfvsafnr(ind),statsKS.ktestrealrembfvsafnr(ind)] = FCDhistktest(FCDvaluesBefore,data.no_remission.afterRest.histFCD,i);
        % real before responders and before non-responders
        [statsKS.hrealrembfvsbfnr(ind),statsKS.prealrembfvsbfnr(ind),statsKS.ktestrealrembfvsbfnr(ind)] = FCDhistktest(data.remission.beforeRest.histFCD,data.no_remission.beforeRest.histFCD,i);
        % real before responders and after responders
        [statsKS.hrealbfrvsafr(ind),statsKS.prealbfrvsafr(ind),statsKS.ktestrealbfrvsafr(ind)] = FCDhistktest(data.remission.afterRest.histFCD,data.remission.beforeRest.histFCD,i);

    end
end

