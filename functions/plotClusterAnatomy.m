function [hLabel] = plotClusterAnatomy(data,labelROI,clust,Order)
% plots paired boxplots for a given measure

c_map =  FMRI_colorMap(8);
hLabel = figure;opt.FontSize = 8;
for cl = 1:clust
        subplot(1,clust,cl)

    for roi = 1:90
	% Plot one single bar as a separate bar series.
	handlebar(roi) = bar(roi, data(cl,Order(roi)), 'BarWidth', 0.6);
	% Apply the color to this bar series.
    if data(cl,Order(roi)) > 0
        set(handlebar(roi), 'FaceColor', [0.6350, 0.0780, 0.1840],'EdgeColor','none');
    elseif data(cl,Order(roi)) < 0
        set(handlebar(roi), 'FaceColor', [0, 0.4470, 0.7410],'EdgeColor','none');
    end
	% Place text atop the bar
	%barTopper = sprintf('y(%d) = %.3f', x(b), eval(b));
	%text(xZscore(b)-0.2, evalZscore(b)+3, char(labelFS{ixZscore}), 'FontSize', barFontSize,'Color',hColor(ixZscore(b),:));
	%set(gca, 'XTickLabels', labelFS{b});
    hold on;
    
    labelROIorder(roi,:) = labelROI(Order(roi),:); % reordering labels
    %title(strcat('Zscore Baseline:', nameCondition,nameGroup))
    
    %if min(evalZscore) > 0
    %else
    %set(gca, 'ylim',[-inf ceil(max(evalZscore))],'XTick',1:7,'XTickLabels', labelFSorderZscore,'FontSize',32);
    %end
    end
    set(gca,'XTick',1:90,'XTickLabels', labelROIorder,'FontSize',8);

    view([90 90])
    %title(measureTitle)
end
 % suptitle('Fractional Occupancy')
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
opt.XMinorTick = 'off';
%opt.YMinorTick = 'off';
hLabel = fancy_figure(hLabel, opt);
end

