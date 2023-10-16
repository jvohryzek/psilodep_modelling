function [H] = FCDdisthist(data,SaveFile,figure_save,fileName,workDate,nameLegend)
% This function calculates the FCD histogram for the different groups
    
    c_map =  FMRI_colorMap(8);

    H = figure
    opt.FontSize = 32;
    subplot(1,2,1)
    histogram([data.remission.beforeRest.histFCD; data.no_remission.beforeRest.histFCD],1000/8,...
        'Normalization','probability','FaceColor',c_map(1,:),'EdgeColor','none')
    hold on
    histogram([data.remission.afterRest.histFCD; data.no_remission.afterRest.histFCD],1000/8,...
        'Normalization','probability','FaceColor',c_map(2,:),'EdgeColor','none')
    legend('before patients','after patients','Location','southwest');
    xlabel('FCD values');ylabel('Probability');
    %xlabel('Cosine Similarity of Leading Eigenvector (FCD values) ');
    subplot(1,2,2)
    histogram(data.remission.beforeRest.histFCD,1000/8,'Normalization','probability','FaceColor',c_map(1,:),'EdgeColor','none')
    hold on
    histogram(data.no_remission.beforeRest.histFCD,1000/8,'Normalization','probability','FaceColor',c_map(2,:),'EdgeColor','none')
    hold on
    histogram(data.remission.afterRest.histFCD,1000/8,'Normalization','probability','FaceColor',c_map(5,:),'EdgeColor','none')
    hold on
    histogram(data.no_remission.afterRest.histFCD,1000/8,'Normalization','probability','FaceColor',c_map(6,:),'EdgeColor','none')
    set(gca,'FontSize',24);
    legend(nameLegend,'Location','southwest');
    %set(hAdd0,'Units', 'Inches', 'Position', [0, 0, 55, 55], 'PaperUnits', 'Inches', 'PaperSize', [55, 55])
    xlabel('FCD values');ylabel('Probability');
    %xlabel('Cosine Similarity of Leading Eigenvector (FCD values) '); ylabel('Probability');
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    H = fancy_figure(H, opt);
    if figure_save
        saveas(H,strcat(SaveFile,workDate,'/PSILODEPresults',fileName,'.jpg'))
        saveas(H,strcat(SaveFile,workDate,'/PSILODEPresults',fileName,'.svg'))
    end
end

