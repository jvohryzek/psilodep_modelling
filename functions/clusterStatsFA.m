function [Cluster_Param] = clusterStatsFA(Clusters,ind,numSbj,numTask,numTp,numClust,respNum,norespNum,typeCase)
% Calculating individual subjects cluster expressions and then dividing
% them into the groups (responders/non-responders)
% fo stands for FRACTIONAL OCCUPANCY
    transferMatrix = zeros(numSbj*numTask,numClust,numClust); % Creating the transfer matrix
    transferMatrixNorm = zeros(numSbj*numTask,numClust,numClust); % Creating the transfer matrix

    % sorting the cluster from the one with most expression
    for sbj = 1:numSbj*numTask
        % 0. FRACTIONAL OCCUPANCY
        % Clusters.IDX => dim(numTpx(numSbjxnumTask),numClust)
        % OPTION 1
        clustExpression = zeros(1,numClust);
        for c = 1:numClust
            clustExpression(c) = mean(Clusters.IDX((1+(numTp*(sbj-1))):(numTp*sbj),1) == c);
        end
        clustExpressionNorm(:,sbj) = clustExpression/sum(clustExpression);
        
        % OPTION 2 - unsuitable when only one cluster takes over
%         clustExpression = hist(Clusters.IDX((1+(numTp*(sbj-1))):(numTp*sbj),1),numClust);
%         clustExpressionNorm(:,sbj) = clustExpression/sum(clustExpression); % divided by number of timepoints
%         clear clustExpression
        
        % TRANSFER MATRIX
        %index = logical((IndexTaskSbj(:,1) == sb).*(IndexTaskSbj(:,2) == tk));
        index = (1+(numTp*(sbj-1))):(numTp*sbj);
        tmp = Clusters.IDX(index);
        % 1. TRANSFER MATRIX
        for tp = 2:numTp
            transferMatrix(sbj,tmp(tp-1),tmp(tp)) = transferMatrix(sbj,tmp(tp-1),tmp(tp)) + 1;
        end
        transferMatrixNorm(sbj,:,:) = squeeze(transferMatrix(sbj,:,:))./sum(squeeze(transferMatrix(sbj,:,:)),2);
        %%% 2. SWITCHING FREQUENCY
        %%% switchFreq(sbj) = mean(diff(tmp)~=0)/TR; % same measure for all clusters
        
        for cl = 1:numClust
            % 3. CLUSTER PROBABILITY
            clustProb(sbj,cl) = mean(tmp == cl); % Like this already normalised!
            % 4. CLUSTER LIFETIME
            clustDur = [];
            clustLifeTime = [];
            clustDur = (tmp == cl);

            % Detect switches in and out of this state
            onset = find(diff([0; clustDur; 0]') == 1);
            offset = find(diff([0; clustDur; 0]') == -1);

            if ~isempty(onset) && ~isempty(offset)
                clustLifeTime = offset - onset;
            else
                clustLifeTime = 0;
            end

            lifeTime(sbj,cl) = mean(clustLifeTime);
        end
    end
    
    switch typeCase
        case 'cross-sectional'
    Cluster_Param{1}.clustExpressionNormSorted = clustExpressionNorm(ind,:); % sorting from the cluster most expressed
    Cluster_Param{1}.clustExpressionNormSorted2 = clustProb(:,ind)'; % sorting from the cluster most expressed
    Cluster_Param{1}.transferMatrixSorted = transferMatrix(:,ind,ind);
    Cluster_Param{1}.transferMatrixNormSorted = transferMatrixNorm(:,ind,ind);
    % defining the groups in the vector of dim(numSbjxnumTask)
    beforeRespSTART = 1; beforeRespEND = respNum;
    afterRespSTART = respNum+1; afterRespEND = 2*respNum;
    beforeNoRespSTART = 2*respNum+1; beforeNoRespEND = 2*respNum+norespNum;
    afterNoRespSTART = 2*respNum+norespNum+1; afterNoRespEND = 2*respNum+2*norespNum;
    
    % real fo values of the groups => Fractional Occupancy 1 and 2 are the same - I calculated
    % it from a different definition
    Cluster_Param{1}.foRestBefore = Cluster_Param{1}.clustExpressionNormSorted(:,[beforeRespSTART:beforeRespEND beforeNoRespSTART:beforeNoRespEND]);
    Cluster_Param{1}.foRestAfter = Cluster_Param{1}.clustExpressionNormSorted(:,[afterRespSTART:afterRespEND afterNoRespSTART:afterNoRespEND]);    
   
    Cluster_Param{1}.foRestBefore2 = Cluster_Param{1}.clustExpressionNormSorted2(:,[beforeRespSTART:beforeRespEND beforeNoRespSTART:beforeNoRespEND]);
    Cluster_Param{1}.foRestAfter2 = Cluster_Param{1}.clustExpressionNormSorted2(:,[afterRespSTART:afterRespEND afterNoRespSTART:afterNoRespEND]);    
    
    Cluster_Param{1}.tmRestBefore = Cluster_Param{1}.transferMatrixSorted([beforeRespSTART:beforeRespEND beforeNoRespSTART:beforeNoRespEND],:,:);
    Cluster_Param{1}.tmRestAfter = Cluster_Param{1}.transferMatrixSorted([afterRespSTART:afterRespEND afterNoRespSTART:afterNoRespEND],:,:);    
    
    Cluster_Param{1}.tmnormRestBefore = Cluster_Param{1}.transferMatrixNormSorted([beforeRespSTART:beforeRespEND beforeNoRespSTART:beforeNoRespEND],:,:);
    Cluster_Param{1}.tmnormRestAfter = Cluster_Param{1}.transferMatrixNormSorted([afterRespSTART:afterRespEND afterNoRespSTART:afterNoRespEND],:,:);    
       
    % real fo values of the groups
    Cluster_Param{1}.foRestBeforeResp = Cluster_Param{1}.clustExpressionNormSorted(:,beforeRespSTART:beforeRespEND);
    Cluster_Param{1}.foRestAfterResp = Cluster_Param{1}.clustExpressionNormSorted(:,afterRespSTART:afterRespEND);
    Cluster_Param{1}.foRestBeforeNoResp = Cluster_Param{1}.clustExpressionNormSorted(:,beforeNoRespSTART:beforeNoRespEND);
    Cluster_Param{1}.foRestAfterNoResp = Cluster_Param{1}.clustExpressionNormSorted(:,afterNoRespSTART:afterNoRespEND);
    
    % diff. fo between groups
    Cluster_Param{1}.foRestDiffResp = Cluster_Param{1}.foRestAfterResp-Cluster_Param{1}.foRestBeforeResp;
    Cluster_Param{1}.foRestDiffNoResp = Cluster_Param{1}.foRestAfterNoResp-Cluster_Param{1}.foRestBeforeNoResp;;

    % mean fo values of the groups
    Cluster_Param{1}.foRestBeforeRespMean = mean(Cluster_Param{1}.clustExpressionNormSorted(:,beforeRespSTART:beforeRespEND),2);
    Cluster_Param{1}.foRestAfterRespMean = mean(Cluster_Param{1}.clustExpressionNormSorted(:,afterRespSTART:afterRespEND),2);
    Cluster_Param{1}.foRestBeforeNoRespMean = mean(Cluster_Param{1}.clustExpressionNormSorted(:,beforeNoRespSTART:beforeNoRespEND),2);
    Cluster_Param{1}.foRestAfterNoRespMean = mean(Cluster_Param{1}.clustExpressionNormSorted(:,afterNoRespSTART:afterNoRespEND),2);
    
    % std. dev fo values of the groups
    Cluster_Param{1}.foRestBeforeRespStd = std(Cluster_Param{1}.clustExpressionNormSorted(:,beforeRespSTART:beforeRespEND),0,2);
    Cluster_Param{1}.foRestAfterRespStd = std(Cluster_Param{1}.clustExpressionNormSorted(:,afterRespSTART:afterRespEND),0,2);
    Cluster_Param{1}.foRestBeforeNoRespStd = std(Cluster_Param{1}.clustExpressionNormSorted(:,beforeNoRespSTART:beforeNoRespEND),0,2);
    Cluster_Param{1}.foRestAfterNoRespStd = std(Cluster_Param{1}.clustExpressionNormSorted(:,afterNoRespSTART:afterNoRespEND),0,2);
    
    %  tm values of the groups
    Cluster_Param{1}.tmRestBeforeResp = Cluster_Param{1}.transferMatrixSorted(beforeRespSTART:beforeRespEND,:,:);
    Cluster_Param{1}.tmRestAfterResp = Cluster_Param{1}.transferMatrixSorted(afterRespSTART:afterRespEND,:,:);
    Cluster_Param{1}.tmRestBeforeNoResp = Cluster_Param{1}.transferMatrixSorted(beforeNoRespSTART:beforeNoRespEND,:,:);    
    Cluster_Param{1}.tmRestAfterNoResp = Cluster_Param{1}.transferMatrixSorted(afterNoRespSTART:afterNoRespEND,:,:);    

    %  tm values of the groups
    Cluster_Param{1}.tmnormRestBeforeResp = Cluster_Param{1}.transferMatrixNormSorted(beforeRespSTART:beforeRespEND,:,:);
    Cluster_Param{1}.tmnormRestAfterResp = Cluster_Param{1}.transferMatrixNormSorted(afterRespSTART:afterRespEND,:,:);
    Cluster_Param{1}.tmnormRestBeforeNoResp = Cluster_Param{1}.transferMatrixNormSorted(beforeNoRespSTART:beforeNoRespEND,:,:);    
    Cluster_Param{1}.tmnormRestAfterNoResp = Cluster_Param{1}.transferMatrixNormSorted(afterNoRespSTART:afterNoRespEND,:,:);    
        
    end     
        
end

