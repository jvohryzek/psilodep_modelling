function [gmeas] = plotClusterFOBars(Cluster_Param)
% plots bars of fractional occupancy for chosen clusters
gmeas = figure; opt.FontSize = 18; subplot(3,3,1)
bar(Cluster_Param{3}{1}.foRestBeforeNoRespMean,'FaceColor',[0 .5 .5]);ylim([0 0.75]);hold on
errorbar([1 2 3],Cluster_Param{3}{1}.foRestBeforeNoRespMean,Cluster_Param{3}{1}.foRestBeforeNoRespStd,...
    '.','LineWidth',3,'Color',[0 .5 .5])
title('Non-Responders before treatment')
set(gca,'xticklabels',{'State1','State2','State3'})
subplot(3,3,3)
bar(Cluster_Param{3}{1}.foRestBeforeRespMean,'FaceColor',[0 .5 .5]);ylim([0 0.75]);hold on
errorbar([1 2 3],Cluster_Param{3}{1}.foRestBeforeRespMean,Cluster_Param{3}{1}.foRestBeforeRespStd,...
    '.','LineWidth',3,'Color',[0 .5 .5])
title('Responders before treatment')
set(gca,'xticklabels',{'State1','State2','State3'})

subplot(3,3,5)
bar(mean(Cluster_Param{3}{1}.foRestBefore,2),'FaceColor',[0 .5 .5]);ylim([0 0.75]);hold on
errorbar([1 2 3],mean(Cluster_Param{3}{1}.foRestBefore,2),std(Cluster_Param{3}{1}.foRestBefore,0,2),...
    '.','LineWidth',3,'Color',[0 .5 .5])
title('All Patients before treatment')
set(gca,'xticklabels',{'State1','State2','State3'})

subplot(3,3,7)
bar(Cluster_Param{3}{1}.foRestAfterNoRespMean,'FaceColor',[0 .5 .5]);ylim([0 0.75]);hold on
errorbar([1 2 3],Cluster_Param{3}{1}.foRestAfterNoRespMean,Cluster_Param{3}{1}.foRestAfterNoRespStd,...
    '.','LineWidth',3,'Color',[0 .5 .5])
title('Non-Responders after treatment')
set(gca,'xticklabels',{'State1','State2','State3'})

subplot(3,3,9)
bar(Cluster_Param{3}{1}.foRestAfterRespMean,'FaceColor',[0 .5 .5]);ylim([0 0.75]);hold on
errorbar([1 2 3],Cluster_Param{3}{1}.foRestAfterRespMean,Cluster_Param{3}{1}.foRestAfterRespStd,...
    '.','LineWidth',3,'Color',[0 .5 .5])
title('Responders after treatment')
set(gca,'xticklabels',{'State1','State2','State3'})

set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
gmeas = fancy_figure(gmeas, opt);

end

