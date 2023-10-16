function [hbrain] = plotClusters2Brain(Vsorted,numClust)
%This function plots all the brain clusters fro the varying cluster value
%sorted based on the cluster's largest fractional occupancy
%
    hbrain = figure;opt.FontSize = 18;

    % Pannel A
    % Plot the cluster centroids over the cortex 
    for i = 3:10
        for m = 1:length(Vsorted{1,i}(:,1))
            subplot(numClust-2,numClust,(i-3)*numClust+m)
            plot_nodes_in_cortex(Vsorted{1,i}(m,:)) % cluster centroids from the largest
            title(['#' num2str(m)])
            axis off
        end
    end
%     subplot(2,5,1);imagesc(tril(hRest,0));title('After/Before'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24)
%     subplot(2,5,2);imagesc(tril(hRestNoResp,0));title('After/Before NoResp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24)
%     subplot(2,5,3);imagesc(tril(hRestResp,0));title('After/Before Resp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24)
%     subplot(2,5,4);imagesc(tril(hRestdiff,0));title(' Resp(After - Before)/NoResp(After - Before)'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24)
%     subplot(2,5,5);imagesc(tril(hRestAfRestNoRest,0));title(' After Resp/NoResp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24) 
%     subplot(2,5,6);imagesc(tril(pRest,0));title('After - Before'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24)
%     subplot(2,5,7);imagesc(tril(pRestNoResp,0));title('After - Before NoResp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24)
%     subplot(2,5,8);imagesc(tril(pRestResp,0));title('After - Before Resp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24)
%     subplot(2,5,9);imagesc(tril(pRestdiff,0));title(' Resp(After - Before)/NoResp(After - Before)'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24) 
%     subplot(2,5,10);imagesc(tril(pRestAfRestNoRest,0));title(' After Resp/NoResp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24) 

    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    hbrain = fancy_figure(hbrain, opt);
end

