function [corr_staticFCVVT, corr_staticFCDFC] = plotAALClusters( Vsorted, numClusters,BOLD_processed_idx,numROI,clustExpressionOrg,Clusters,ind,staticFC_mean, Order )
colormap(jet)    
VVT_mean = zeros(numROI);
dFC_mean = zeros(numROI,numROI);
for cl = 1:numClusters
    % Pannel A
    % Plot the cluster centroids over the cortex 
    subplot(3,numClusters+1,cl)
    plot_nodes_in_cortex(Vsorted(cl,:)) % cluster centroids from the largest
    title(['#' num2str(cl)])
    axis off
    % Pannel C 
    % Plot the centroids outerproduct
    subplot(3,numClusters+1,cl+numClusters+1)
    VVT = Vsorted(cl,:)'*Vsorted(cl,:);
    colormap(jet)

    imagesc(VVT(Order,Order),[-max(max(abs(VVT))),max(max(abs(VVT)))])   % the order has got to do with the brain regions?
    axis square
    title('Outer product') 
    ylabel('Brain area #')
    xlabel('Brain area #') 
    colorbar;
    VVT_mean = VVT_mean + VVT*clustExpressionOrg(cl);
    % Panel E
    % Plot dFC states averaged (static dFC states)
    subplot(3,numClusters+1,cl+2*numClusters+2)
    dFC = corrcoef(BOLD_processed_idx(:,Clusters.IDX == ind(cl))'); % making sure clusters are in order
    imagesc(dFC(Order,Order))   % the order has got to do with the brain regions?
    axis square
    title('dFC states') 
    ylabel('Brain area #')
    xlabel('Brain area #') 
    colorbar;% caxis([-0.5 0.05]);
    dFC_mean = dFC_mean + dFC*clustExpressionOrg(cl);
    dFC_mean_selen(:,:,cl) = dFC_mean;

end
VVT_mean = VVT_mean/sum(clustExpressionOrg);
dFC_mean = dFC_mean/sum(clustExpressionOrg);

% Pannel B
subplot(3,numClusters+1,numClusters+1)
imagesc(staticFC_mean(Order,Order))
title({'Static FC'})
ylabel('Brain area')
xlabel('Brain area')
colorbar;% caxis([-0.5 0.05]);
axis square
% Pannel D plot the weighted sum of outer products

subplot(3,numClusters+1,(numClusters+1)*2)
imagesc(VVT_mean(Order,Order))
title({'Weighted sum of VV^T'})
ylabel('Brain area')
xlabel('Brain area')
colorbar
axis square
% Panel F dFC_mean
subplot(3,numClusters+1,(numClusters+1)*3)
imagesc(dFC_mean(Order,Order))
title({'Weighted sum of dFC'})
ylabel('Brain area')
xlabel('Brain area')
colorbar;% caxis([-0.5 0.05]);
axis square
corr_staticFCVVT = corrcoef(VVT_mean,staticFC_mean);
corr_staticFCDFC = corrcoef(dFC_mean,staticFC_mean);
end

