function [ cl ] = plotClusters( V, numClusters,BOLD_processed_idx,numROI,clustExpressionOrg,Clusters,ind,staticFC_mean )
%UNTITLED10 Summary of this function goes here
colormap(jet)    
VVT_mean = zeros(numROI);
dFC_mean = zeros(numROI,numROI);
for cl = 1:numClusters
    % Pannel A
    % Plot the cluster centroids over the cortex 
    subplot(3,numClusters+1,cl)
    NodeMatrixPlots(V(cl,:)) % cluster centroids from the largest
    title(['#' num2str(cl)])
    axis off
    % Pannel C 
    % Plot the centroids outerproduct
    subplot(3,numClusters+1,cl+numClusters+1)
    VVT = V(cl,:)'*V(cl,:);   
    imagesc(VVT)   % the order has got to do with the brain regions?
    axis square
    title('Outer product') 
    ylabel('Brain area #')
    xlabel('Brain area #') 
    VVT_mean = VVT_mean + VVT*clustExpressionOrg(cl);
    % Panel E
    % Plot dFC states averaged (static dFC states)
    subplot(3,numClusters+1,cl+2*numClusters+2)
    dFC = corrcoef(BOLD_processed_idx(:,Clusters.IDX == ind(cl))'); % making sure clusters are in order
    imagesc(dFC)   % the order has got to do with the brain regions?
    axis square
    title('dFC states') 
    ylabel('Brain area #')
    xlabel('Brain area #') 
    dFC_mean = dFC_mean + dFC*clustExpressionOrg(cl);
    dFC_mean_selen(:,:,cl) = dFC_mean;

end
VVT_mean = VVT_mean/sum(clustExpressionOrg);
dFC_mean = dFC_mean/sum(clustExpressionOrg);

% Pannel B
subplot(3,numClusters+1,numClusters+1)
imagesc(staticFC_mean)
title({'Static FC'})
ylabel('Brain area')
xlabel('Brain area')
colorbar
axis square
% Pannel D plot the weighted sum of outer products

subplot(3,numClusters+1,(numClusters+1)*2)
imagesc(VVT_mean)
title({'Weighted sum of VV^T'})
ylabel('Brain area')
xlabel('Brain area')
colorbar
axis square
% Panel F dFC_mean
subplot(3,numClusters+1,(numClusters+1)*3)
imagesc(dFC_mean)
title({'Weighted sum of dFC'})
ylabel('Brain area')
xlabel('Brain area')
colorbar
axis square


end

