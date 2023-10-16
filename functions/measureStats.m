function [Measure_Param,Measure_Stats] = measureStats(data,measureType,ConfInt,tailInt,numPerm,indTMm,indTMn)
% Calculating individual subjects measures statistics (responders/non-responders)
%
switch measureType
    case 'metastability'
        dataRestBeforeResp = data.remission.beforeRest.Metastability;
        dataRestAfterResp = data.remission.afterRest.Metastability;
        dataRestBeforeNoResp = data.no_remission.beforeRest.Metastability;
        dataRestAfterNoResp = data.no_remission.afterRest.Metastability;
    case 'synchrony'
        dataRestBeforeResp = data.remission.beforeRest.Synchrony;
        dataRestAfterResp = data.remission.afterRest.Synchrony;
        dataRestBeforeNoResp = data.no_remission.beforeRest.Synchrony;
        dataRestAfterNoResp = data.no_remission.afterRest.Synchrony;
    case 'switchFreq'
        dataRestBeforeResp = data.remission.beforeRest.switchFreq;
        dataRestAfterResp = data.remission.afterRest.switchFreq;
        dataRestBeforeNoResp = data.no_remission.beforeRest.switchFreq;
        dataRestAfterNoResp = data.no_remission.afterRest.switchFreq;
   case 'lifeTime'
        dataRestBeforeResp = data.remission.beforeRest.lifeTime(:,indTMm)';
        dataRestAfterResp = data.remission.afterRest.lifeTime(:,indTMm)';
        dataRestBeforeNoResp = data.no_remission.beforeRest.lifeTime(:,indTMm)';
        dataRestAfterNoResp = data.no_remission.afterRest.lifeTime(:,indTMm)';
    case 'transferMatrix'
        dataRestBeforeResp = squeeze(data.remission.beforeRest.transferMatrix(indTMm,indTMn,:))';
        dataRestAfterResp = squeeze(data.remission.afterRest.transferMatrix(indTMm,indTMn,:))';
        dataRestBeforeNoResp = squeeze(data.no_remission.beforeRest.transferMatrix(indTMm,indTMn,:))';
        dataRestAfterNoResp = squeeze(data.no_remission.afterRest.transferMatrix(indTMm,indTMn,:))';
    end
    % real fo values of the groups
    Measure_Param.RestBefore = [dataRestBeforeResp dataRestBeforeNoResp];
    Measure_Param.RestAfter = [dataRestAfterResp dataRestAfterNoResp];    
    
    % real fo values of the groups
    Measure_Param.RestBeforeResp = dataRestBeforeResp;
    Measure_Param.RestAfterResp = dataRestAfterResp;
    Measure_Param.RestBeforeNoResp = dataRestBeforeNoResp;
    Measure_Param.RestAfterNoResp = dataRestAfterNoResp;
    
    % diff. fo between groups
    Measure_Param.RestDiffResp = dataRestAfterResp-dataRestBeforeResp;
    Measure_Param.RestDiffNoResp = dataRestAfterNoResp-dataRestBeforeNoResp;

    % mean fo values of the groups
    Measure_Param.RestBeforeRespMean = mean(dataRestBeforeResp,2);
    Measure_Param.RestAfterRespMean = mean(dataRestAfterResp,2);
    Measure_Param.RestBeforeNoRespMean = mean(dataRestBeforeNoResp,2);
    Measure_Param.RestAfterNoRespMean = mean(dataRestAfterNoResp,2);
    % std. dev fo values of the groups
    Measure_Param.RestBeforeRespStd = std(dataRestBeforeResp,0,2);
    Measure_Param.RestAfterRespStd = std(dataRestAfterResp,0,2);
    Measure_Param.RestBeforeNoRespStd = std(dataRestBeforeNoResp,0,2);
    Measure_Param.RestAfterNoRespStd = std(dataRestAfterNoResp,0,2);    
    
    %% statistics
    
%         % ttest
%         [Measure_Stats.hRest,Measure_Stats.pRest] = ttest(Measure_Param.RestAfter,Measure_Param.RestBefore,'Alpha',ConfInt,'tail',tailInt);
%         [Measure_Stats.hRestResp,Measure_Stats.pRestResp] = ttest(Measure_Param.RestAfterResp,Measure_Param.RestBeforeResp,'Alpha',ConfInt,'tail',tailInt);
%         [Measure_Stats.hRestNoResp,Measure_Stats.pRestNoResp] = ttest(Measure_Param.RestAfterNoResp,Measure_Param.RestBeforeNoResp,'Alpha',ConfInt,'tail',tailInt);
%         [Measure_Stats.hRestdiff,Measure_Stats.pRestdiff] = ttest2(Measure_Param.RestDiffResp,Measure_Param.RestDiffNoResp,'Alpha',ConfInt,'tail',tailInt); % one sided
%         [Measure_Stats.hRestAfRestNoRest,Measure_Stats.pRestAfRestNoRest] = ttest2(Measure_Param.RestAfterResp,Measure_Param.RestAfterNoResp,'Alpha',ConfInt,'tail',tailInt); % one sided
%         [Measure_Stats.hRestBfAllAfRest,Measure_Stats.pRestBfAllAfRest] = ttest2(Measure_Param.RestAfterResp,Measure_Param.RestBefore,'Alpha',ConfInt,'tail',tailInt); % one sided
%         % ranksum        
%         [Measure_Stats.pRestWsignrank,Measure_Stats.hRestWsignrank] = signrank(Measure_Param.RestAfter,Measure_Param.RestBefore,'alpha',ConfInt,'tail',tailInt);
%         [Measure_Stats.pRestRespWsignrank,Measure_Stats.hRestRespWsignrank] = signrank(Measure_Param.RestAfterResp,Measure_Param.RestBeforeResp,'alpha',ConfInt,'tail',tailInt);
%         [Measure_Stats.pRestNoRespWsignrank,Measure_Stats.hRestNoRespWsignrank] = signrank(Measure_Param.RestAfterNoResp,Measure_Param.RestBeforeNoResp,'alpha',ConfInt,'tail',tailInt);
%         [Measure_Stats.pRestdiffWranksum,Measure_Stats.hRestdiffWranksum] = ranksum(Measure_Param.RestDiffResp,Measure_Param.RestDiffNoResp,'alpha',ConfInt,'tail','left'); % one sided
%         [Measure_Stats.pRestAfRestNoRestWranksum,Measure_Stats.hRestAfRestNoRestWranksum] = ranksum(Measure_Param.RestAfterResp,Measure_Param.RestAfterNoResp,'alpha',ConfInt,'tail',tailInt); % one sided
%         [Measure_Stats.pRestBfAllAfRestWranksum,Measure_Stats.hRestBfAllAfRestWranksum] = ranksum(Measure_Param.RestAfterResp,Measure_Param.RestBefore,'alpha',ConfInt,'tail',tailInt); % one sided
        % permuted ttest
        disp('Compare Cluster Probability of before All with after All')
        statsRest = permutation_htesting_np([Measure_Param.RestBefore,Measure_Param.RestAfter],...
            [ones(1,numel(Measure_Param.RestBefore)) 2*ones(1,numel(Measure_Param.RestAfter))],numPerm,0.05,'ttest');
        Measure_Stats.pRestPermTtest = min(statsRest.pvals);
        Measure_Stats.hRestPermTtest = min(statsRest.pvals) < 0.05;

        disp('Compare Cluster Probability of before Resp with after Resp')
        statsRestResp = permutation_htesting_np([Measure_Param.RestAfterResp,Measure_Param.RestBeforeResp],...
            [ones(1,numel(Measure_Param.RestAfterResp)) 2*ones(1,numel(Measure_Param.RestBeforeResp))],numPerm,0.05,'ttest');
        Measure_Stats.pRestRespPermTtest = min(statsRestResp.pvals);
        Measure_Stats.hRestRespPermTtest = min(statsRestResp.pvals) < 0.05;

        disp('Compare Cluster Probability of before NoResp with after NoResp')
        statsRestNoResp = permutation_htesting_np([Measure_Param.RestAfterNoResp,Measure_Param.RestBeforeNoResp],...
            [ones(1,numel(Measure_Param.RestAfterNoResp)) 2*ones(1,numel(Measure_Param.RestBeforeNoResp))],numPerm,0.05,'ttest');
        Measure_Stats.pRestNoRespPermTtest = min(statsRestNoResp.pvals);
        Measure_Stats.hRestNoRespPermTtest = min(statsRestNoResp.pvals) < 0.05;

        disp('Compare Cluster Probability of Resp Diff with NoResp Diff')
        statsRestdiff = permutation_htesting_np([Measure_Param.RestDiffResp,Measure_Param.RestDiffNoResp],...
            [ones(1,numel(Measure_Param.RestDiffResp)) 2*ones(1,numel(Measure_Param.RestDiffNoResp))],numPerm,0.05,'ttest2');
        Measure_Stats.pRestdiffPermTtest2 = min(statsRestdiff.pvals);
        Measure_Stats.hRestdiffPermTtest2 = min(statsRestdiff.pvals) < 0.05;

        disp('Compare Cluster Probability of before Resp with after Resp')
        statsRestAfRestNoRest = permutation_htesting_np([Measure_Param.RestAfterResp,Measure_Param.RestAfterNoResp],...
            [ones(1,numel(Measure_Param.RestAfterResp)) 2*ones(1,numel(Measure_Param.RestAfterNoResp))],numPerm,0.05,'ttest2');
        Measure_Stats.pRestAfRestNoRestPermTtest2 = min(statsRestAfRestNoRest.pvals);
        Measure_Stats.hRestAfRestNoRestPermTtest2 = min(statsRestAfRestNoRest.pvals) < 0.05;

        disp('Compare Cluster Probability of After NoResp with before All')
        statsRestBfAllAfRest = permutation_htesting_np([Measure_Param.RestAfterResp,Measure_Param.RestBefore],...
            [ones(1,numel(Measure_Param.RestAfterResp)) 2*ones(1,numel(Measure_Param.RestBefore))],numPerm,0.05,'ttest2');
        Measure_Stats.pRestBfAllAfRestPermTtest2 = min(statsRestBfAllAfRest.pvals);
        Measure_Stats.hRestBfAllAfRestPermTtest2 = min(statsRestBfAllAfRest.pvals) < 0.05;

        % permuted ranksum

        disp('Compare Cluster Probability of before All with after All')
        statsRest = permutation_htesting_np([Measure_Param.RestBefore,Measure_Param.RestAfter],...
            [ones(1,numel(Measure_Param.RestBefore)) 2*ones(1,numel(Measure_Param.RestAfter))],numPerm,0.05,'signrank');
        Measure_Stats.pRestPermSignrank = min(statsRest.pvals);
        Measure_Stats.hRestPermSignrank = min(statsRest.pvals) < 0.05;

        disp('Compare Cluster Probability of before Resp with after Resp')
        statsRestResp = permutation_htesting_np([Measure_Param.RestAfterResp,Measure_Param.RestBeforeResp],...
            [ones(1,numel(Measure_Param.RestAfterResp)) 2*ones(1,numel(Measure_Param.RestBeforeResp))],numPerm,0.05,'signrank');
        Measure_Stats.pRestRespPermSignrank = min(statsRestResp.pvals);
        Measure_Stats.hRestRespPermSignrank = min(statsRestResp.pvals) < 0.05;

        disp('Compare Cluster Probability of before NoResp with after NoResp')
        statsRestNoResp = permutation_htesting_np([Measure_Param.RestAfterNoResp,Measure_Param.RestBeforeNoResp],...
            [ones(1,numel(Measure_Param.RestAfterNoResp)) 2*ones(1,numel(Measure_Param.RestBeforeNoResp))],numPerm,0.05,'signrank');
        Measure_Stats.pRestNoRespPermSignrank = min(statsRestNoResp.pvals);
        Measure_Stats.hRestNoRespPermSignrank = min(statsRestNoResp.pvals) < 0.05;

        disp('Compare Cluster Probability of Resp Diff with NoResp Diff')
        statsRestdiff = permutation_htesting_np([Measure_Param.RestDiffResp,Measure_Param.RestDiffNoResp],...
            [ones(1,numel(Measure_Param.RestDiffResp)) 2*ones(1,numel(Measure_Param.RestDiffNoResp))],numPerm,0.05,'ranksum');
        Measure_Stats.pRestdiffPermRanksum = min(statsRestdiff.pvals);
        Measure_Stats.hRestdiffPermRanksum = min(statsRestdiff.pvals) < 0.05;

        disp('Compare Cluster Probability of before Resp with after Resp')
        statsRestAfRestNoRest = permutation_htesting_np([Measure_Param.RestAfterResp,Measure_Param.RestAfterNoResp],...
            [ones(1,numel(Measure_Param.RestAfterResp)) 2*ones(1,numel(Measure_Param.RestAfterNoResp))],numPerm,0.05,'ranksum');
        Measure_Stats.pRestAfRestNoRestPermRanksum = min(statsRestAfRestNoRest.pvals);
        Measure_Stats.hRestAfRestNoRestPermRanksum = min(statsRestAfRestNoRest.pvals) < 0.05;

        disp('Compare Cluster Probability of After NoResp with before All')
        statsRestBfAllAfRest = permutation_htesting_np([Measure_Param.RestAfterResp,Measure_Param.RestBefore],...
            [ones(1,numel(Measure_Param.RestAfterResp)) 2*ones(1,numel(Measure_Param.RestBefore))],numPerm,0.05,'ranksum');
        Measure_Stats.pRestBfAllAfRestPermRanksum = min(statsRestBfAllAfRest.pvals);
        Measure_Stats.hRestBfAllAfRestPermRanksum = min(statsRestBfAllAfRest.pvals) < 0.05;
    
    
end

