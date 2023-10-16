function [ staticFC_Yeo, staticFC_Yeo_with_conn, staticFC_Yeo_betw_conn] = sFC4YeoAnalysis( staticFC,ROI2Yeo)
% calculating within and between network connectivity
labelYeo = {'VisCent','VisPeri','SomMotA','SomMotB','DorsAttnA','DorsAttnB',...
    'SalVentAttnA','SalVentAttnB','LimbicA','LimbicB','ContA','ContB','ContC',...
    'DefaultA','DefaultB','DefaultC','TempPar'};
Isubdiag_Yeo = find(tril(ones(17),-1));

for roi=1:17
    idxYeo(roi,:) = contains(ROI2Yeo,labelYeo{roi});
end
for i=1:17
    for j=1:17
        staticFC_Yeo(i,j) = mean(mean(staticFC(idxYeo(i,:),idxYeo(j,:))));
        
    end
end 
staticFC_Yeo_with_conn = diag(staticFC_Yeo);
staticFC_Yeo_betw_conn = staticFC_Yeo(Isubdiag_Yeo);
end

