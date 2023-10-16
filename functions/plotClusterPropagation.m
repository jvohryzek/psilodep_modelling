function [hPropagation] = plotClusterPropagation(data,labelROI,clust,Order)
% plots paired boxplots for a given measure

c_map =  FMRI_colorMap(8);
hPropagation = figure;opt.FontSize = 8;
data = data';

for cl = 1:length(data(:,1))
        subplot(length(data(:,1))/(length(data(:,1))/3),(length(data(:,1))/3),cl)
    histData = histogram(max( data));
    histEd = histData.BinEdges;
    for roi = 1:90
	% Plot one single bar as a separate bar series.
	handlebar(roi) = bar(roi, data(cl,Order(roi)), 'BarWidth', 0.6);
	% Apply the color to this bar series.
    if histEd(end-1) < data(cl,Order(roi))
        set(handlebar(roi), 'FaceColor', [0.6350, 0.0780, 0.1840],'EdgeColor','none','FaceAlpha', 1);
    elseif histEd(end-2) < data(cl,Order(roi)) && data(cl,Order(roi)) < histEd(end-1)
        set(handlebar(roi), 'FaceColor', [0.6350, 0.0780, 0.1840],'EdgeColor','none','FaceAlpha', 1);
        %alpha(0.8) [0.6350, 0.2780, 0.1840]
    elseif histEd(end-3) < data(cl,Order(roi)) && data(cl,Order(roi)) < histEd(end-2)
        set(handlebar(roi), 'FaceColor', [0.6350, 0.0780, 0.1840],'EdgeColor','none','FaceAlpha', 1);
        %alpha(0.6) [0.6350, 0.4780, 0.1840]
    elseif 0 < data(cl,Order(roi)) && data(cl,Order(roi)) < histEd(end-3)
        set(handlebar(roi), 'FaceColor', [0.6350, 0.0780, 0.1840],'EdgeColor','none','FaceAlpha', 1);
        %alpha(0.4) [0.6350, 0.6780, 0.1840]
    elseif data(cl,Order(roi)) < 0
        set(handlebar(roi), 'FaceColor', [0, 0.4470, 0.7410],'EdgeColor','none');
    end

	% Place text atop the bar
	%barTopper = sprintf('y(%d) = %.3f', x(b), eval(b));
	%text(xZscore(b)-0.2, evalZscore(b)+3, char(labelFS{ixZscore}), 'FontSize', barFontSize,'Color',hColor(ixZscore(b),:));
	%set(gca, 'XTickLabels', labelFS{b});
     hold on;
    
    % labelROIorder(roi,:) = labelROI(Order(roi),:); % reordering labels
    %title(strcat('Zscore Baseline:', nameCondition,nameGroup))
    
    %if min(evalZscore) > 0
    %else
    %set(gca, 'ylim',[-inf ceil(max(evalZscore))],'XTick',1:7,'XTickLabels', labelFSorderZscore,'FontSize',32);
    %end
    end
    % set(gca,'XTick',1:90,'XTickLabels', labelROIorder,'FontSize',8);
    axis off
    view([90 90])
    %title(measureTitle)
end
 % suptitle('Fractional Occupancy')
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
opt.XMinorTick = 'off';
%opt.YMinorTick = 'off';
hPropagation = fancy_figure(hPropagation, opt);
end

