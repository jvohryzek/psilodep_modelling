function [h3,corr_staticFCVVT, corr_staticFCDFC] = plotClusters2BrainStats( Vsorted, numClust,BOLD_processed_idx,numAreas,clustExpressionOrg,Clusters,ind,staticFC_mean,Order)
% plots the cluster on a brain, static FC, elading eigenvectors and more
% consult plotAALClusters
%
        h3 = figure;opt.FontSize = 18;
        [corr_staticFCVVT, corr_staticFCDFC] = plotAALClusters(Vsorted, numClust,BOLD_processed_idx,numAreas,clustExpressionOrg,Clusters,ind,staticFC_mean,Order )
        set(gcf, 'units','normalized','outerposition',[0 0 2 2]);
        h3 = fancy_figure(h3, opt);
end

