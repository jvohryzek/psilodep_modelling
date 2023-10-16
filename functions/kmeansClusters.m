function [Kmeans_results, distM_FCD, ind_max, dunnScore, evaDB, evaSilh] = kmeansClusters(Leading_Eig_idx, maxK, replicates, distance )
% Kmeans clustering
% opt= statset('UseParallel',1); %,'UseSubstreams',1);
% The options may vary according to the Matlab version
Kmeans_results = cell(1,maxK);
for k = 2:maxK
    
    disp(['Calculating for ' num2str(k) 'clusters'])
    [IDX, C, SUMD, D] = kmeans(Leading_Eig_idx',k,'Distance',distance,'Replicates',replicates,'Display','final'); %,'Options',opt);   
    clust(:,k-1) = IDX;
    Kmeans_results{k}.IDX = IDX;
    Kmeans_results{k}.C = C; 
    Kmeans_results{k}.SUMD = SUMD; 
    Kmeans_results{k}.D = D;
end

%% Evaluate Clustering performance
% 1. Dunn's score
if strcmp(distance,'sqeuclidean')
   distance = 'euclidean';
elseif strcmp(distance,'cityblock')
   distance = 'cityblock';
elseif strcmp(distance,'correlation')
   distance = 'correlation';
elseif strcmp(distance,'cosine')
   distance = 'cosine';
end
distM_FCD = squareform(pdist(Leading_Eig_idx',distance));
dunnScore = zeros(maxK,1);
for j=2:maxK
    dunnScore(j)=dunns(j,distM_FCD,Kmeans_results{j}.IDX);
    disp(['Performance for ' num2str(j) 'clusters'])
end
[~,ind_max] = max(dunnScore);
disp(['Best clustering solution: ' num2str(ind_max) ' clusters']);

% 2. Silhouette

if strcmp(distance,'euclidean')
    eva = evalclusters(Leading_Eig_idx',clust,'Silhouette','Distance','Euclidean')
    evaSilh = eva.CriterionValues;
elseif strcmp(distance,'cityblock')
   eva = evalclusters(Leading_Eig_idx',clust,'Silhouette','Distance','cityblock');
   evaSilh = eva.CriterionValues;
elseif strcmp(distance,'cosine')
   eva = evalclusters(Leading_Eig_idx',clust,'Silhouette','Distance','cosine');
   evaSilh = eva.CriterionValues;
elseif strcmp(distance,'correlation')
   eva = evalclusters(Leading_Eig_idx',clust,'Silhouette','Distance','correlation');
   evaSilh = eva.CriterionValues;
end

% 3. Davies-Bouldin => for euclidean distance only
% [1] Davies, D. L., and D. W. Bouldin. ?A Cluster Separation Measure.? IEEE Transactions on Pattern Analysis and Machine Intelligence. Vol. PAMI-1, No. 2, 1979, pp. 224?227.

if strcmp(distance,'euclidean')

    eva = evalclusters(Leading_Eig_idx',clust,'DaviesBouldin');
    evaDB = eva.CriterionValues;
% elseif strcmp(distance,'cityblock')
% 
%    eva = evalclusters(Leading_Eig_idx',clust,'DaviesBouldin');
%    evaDB = eva.CriterionValues;
else
    evaDB = [];
end

end
