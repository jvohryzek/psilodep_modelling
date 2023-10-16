function [h2] = kmeansClustersEval(distance,dunnScore,evaDB,evaSilh)
% Evaluates kemans cluster parameters
h2 = figure;opt.FontSize = 18;
if strcmp(distance,'sqeuclidean')
subplot(1,3,1); plot(2:10,dunnScore(2:end));
title('Dunns Score'),xlabel('# of clusters'),ylabel('Dunns Score')
subplot(1,3,2); plot(2:10,evaDB)
title('Davies-Bouldin Score'),xlabel('# of clusters'),ylabel('Davies-Bouldin Score')
subplot(1,3,3); plot(2:10,evaSilh)
title('Silhouette Score'),xlabel('# of clusters'),ylabel('Silhouette Score') 

elseif strcmp(distance,'cityblock')
subplot(1,2,1); plot(2:10,dunnScore(2:end));
title('Dunns Score'),xlabel('# of clusters'),ylabel('Dunns Score')
subplot(1,2,2); plot(2:10,evaSilh)
title('Silhouette Score'),xlabel('# of clusters'),ylabel('Silhouette Score') 

elseif strcmp(distance,'cosine')
subplot(1,2,1); plot(2:10,dunnScore(2:end));
title('Dunns Score'),xlabel('# of clusters'),ylabel('Dunns Score')
subplot(1,2,2); plot(2:10,evaSilh)
title('Silhouette Score'),xlabel('# of clusters'),ylabel('Silhouette Score') 

elseif strcmp(distance,'correlation')
subplot(1,2,1); plot(2:10,dunnScore(2:end));
title('Dunns Score'),xlabel('# of clusters'),ylabel('Dunns Score')
subplot(1,2,2); plot(2:10,evaSilh)
title('Silhouette Score'),xlabel('# of clusters'),ylabel('Silhouette Score') 
end
suptitle(distance)
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
h2 = fancy_figure(h2, opt);
end

