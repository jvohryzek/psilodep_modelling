function [hmeasBar] = plotMeasureBars(data,dataStats,numClust,measureTitle)
% plots paired boxplots for a given measure

c_map =  FMRI_colorMap(8);
hmeasBar = figure;opt.FontSize = 18;
for jing = 1:2
    subplot(1,2,jing)
    if jing == 1
        tmpY = data.RestBeforeResp';
        tmpO = data.RestAfterResp';

    elseif jing == 2
         tmpY = data.RestBeforeNoResp';
         tmpO = data.RestAfterNoResp';

    end
    boxplot([tmpY tmpO])%,'outliersize' ,0.001);%, [ones(size(tmp_SD_last)) ones(size(tmp_SD_first))*2]);
    h = findobj(gca,'Tag','Box');
    R1 = [0,150/255,100/255];
    R2 = c_map(2,:); %[44/255,139/255,116/255];
    R3 = [0 0.5 0]; %[230/255,160/255,70/255];

    for j=1:2
       patch(get(h(j),'XData'),get(h(j),'YData'),R1,'FaceAlpha',.4, 'LineWidth', 2, 'EdgeColor', R1,'EdgeAlpha',0.7);
    end
    patch(get(h(1),'XData'),get(h(1),'YData'),R3,'FaceAlpha',.4, 'LineWidth', 2, 'EdgeColor', R3,'EdgeAlpha',0.7);
    set(gca, 'Xtick', [1:2], 'Xticklabel', {'before','after'},'FontSize',20);

    R_asc = [0.6,0.6,0.6];
    R_des = [0.8,0.8,0.8];

    for j=1:length(tmpY)
        if tmpY(j)>tmpO(j); R = R_asc; else R = R_des; end;
        hold on;
        rand_dist = 0.1*(rand-0.5);
        plot([1+rand_dist 2+rand_dist], [tmpY(j) tmpO(j)], '-.', 'Color', R, 'MarkerSize', 1,'LineWidth',2);
        hold on;
        plot([1+rand_dist 2+rand_dist], [tmpY(j) tmpO(j)], '.', 'MarkerEdgeColor', R, 'MarkerSize', 20);
    end
    
    if jing ==1
        ylabel(strcat(measureTitle,' Resp '))
        title(strcat(measureTitle,': p = ',num2str(dataStats.pRestRespPermTtest)))
    elseif jing == 2
        ylabel(strcat(measureTitle,' Non-Resp '))
        title(strcat(measureTitle,': p = ',num2str(dataStats.pRestNoRespPermTtest)))
    end
    
end
 % suptitle('Fractional Occupancy')
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
hmeasBar = fancy_figure(hmeasBar, opt);
%title(measureTitle)
end

