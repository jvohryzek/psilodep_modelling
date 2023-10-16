function [clustExpressionOrg,ind,Vsorted] = kmeansClustersSorted(Clusters,numClust,typeSort)
% SORTING CLUSTER BASED ON THEIR EXPRESSION IN TIME
    clustExpression           = hist(Clusters.IDX,numClust);
    [clustExpressionOrg, ind] = sort(clustExpression,typeSort);
    Vsorted                   = Clusters.C(ind,:); % sorted Cluster Centroids
end

