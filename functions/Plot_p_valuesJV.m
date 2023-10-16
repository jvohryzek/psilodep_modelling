function [F1, F2, F3] = Plot_p_valuesJV(P_pval,rangeK,Kmeans_results)
%load LEiDA_results.mat P_pval LT_pval rangeK Kmeans_results
% 25/03/2019: adjusted in relation to the start brain plots but unlike lausanne no
% brain half so just to keep the code but use the normal function
mink=rangeK(1);
maxk=rangeK(end);

%% PROBABILITIES

sigC=zeros(1,length(rangeK));

for k=3:length(rangeK)   
    [Min_p_value(k), sigC(k)]=min(P_pval(k,P_pval(k,:)>0));    
end

[F1] = figure
semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
hold on
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'g--','LineWidth',1)
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./sum(rangeK)*ones(1,length(rangeK)+2),'b--','LineWidth',1)

for k=1:length(rangeK) 
    for c=1:rangeK(k)        
        if P_pval(k,c)>0.05
            semilogy(rangeK(k),P_pval(k,c),'*k');
        end
        if P_pval(k,c)<0.05 && P_pval(k,c)>(0.05/rangeK(k))
            semilogy(rangeK(k),P_pval(k,c),'*r');
        end
        if P_pval(k,c)<(0.05/rangeK(k)) && P_pval(k,c)>(0.05/sum(rangeK))
            semilogy(rangeK(k),P_pval(k,c),'*g');
        end
        if P_pval(k,c)<=(0.05/sum(rangeK))
            semilogy(rangeK(k),P_pval(k,c),'*b');
        end
    end
end

ylabel('Prob Patients vs Controls (p-value)')
xlabel('Number of clusters K')
xlim([rangeK(1)-1 rangeK(end)+1])
box off

[F2] = figure('Name','Probabilities')
load AAL_labels.mat label90
N_areas=90;
Order=[1:2:N_areas N_areas:-2:2];
label90=label90(Order,:);

V_all=zeros(N_areas,length(rangeK));

for k=3:length(rangeK)
    c=sigC(k);
    V=Kmeans_results{k}(c,:) % Kmeans_results{k}.C(c,:);
    V=V/max(abs(V));
    V_all(:,k)=V;
    
    V=V(Order);
    subplot(1,length(rangeK),k)
    hold on
    barh((V.*(V<=0)),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
    barh((V.*(V>0)),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
    ylim([0 N_areas+1])
    xlim([-1 1])
    
    grid on
    set(gca,'YTick',1:N_areas,'Fontsize',8)    
    if k==1
        set(gca,'YTickLabel',label90(end:-1:1,:))
    else
        set(gca,'YTickLabel',[])
    end
    if Min_p_value(k)<0.05 && Min_p_value(k)>(0.05/rangeK(k))
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k))]},'Fontsize',10,'color','r')
    elseif Min_p_value(k)>0.05
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k))]},'Fontsize',10,'color','k')
    elseif Min_p_value(k)<(0.05/k) && Min_p_value(k)>(0.05/sum(rangeK))
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k))]},'Fontsize',10,'color','g')
    elseif Min_p_value(k)<(0.05/sum(rangeK))
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k))]},'Fontsize',10,'color','b')        
    end
end

V_best=Kmeans_results{rangeK==6}(3,:)
% V_best=Kmeans_results{rangeK==6}.C(3,:);

[F3] = figure
subplot(3,4,1)
plot_nodes_in_cortex(V_best)
rotate3d
view(0,0)
title('K=6, c=3')

subplot(3,4,3)
plot_nodes_in_cortex(V_best)
view(-90,0)
rotate3d
title('K=6, c=3')
% pushing the plot
pos = get(gca, 'Position');
pos(1) = pos(1)-0.1;
set(gca, 'Position', pos)

subplot(3,4,6)
plot_nodes_in_cortex(V_best)
view(-90,90)
rotate3d
% pushing the plot
pos = get(gca, 'Position');
pos(1) = pos(1)-0.05;
pos(2) = pos(2)+0.05;
set(gca, 'Position', pos)

subplot(3,4,9)
plot_nodes_in_cortex(V_best)
view(0,90)
rotate3d

% pushing the plot
pos = get(gca, 'Position');
pos(2) = pos(2)+0.1;
set(gca, 'Position', pos)

subplot(3,4,11)
plot_nodes_in_cortex(V_best)
view(0,90)
rotate3d
% pushing the plot
pos = get(gca, 'Position');
pos(1) = pos(1)-0.1;
pos(2) = pos(2)+0.1;
set(gca, 'Position', pos)

V_best=V_best(Order);
V_best=V_best./max(abs(V_best));

subplot(3,4,[4 8 12])
hold on
barh(V_best.*(V_best<=0),'FaceColor','blue','EdgeColor','none')
barh(V_best.*(V_best>0),'FaceColor','red','EdgeColor','none','Barwidth',.4)
set(gca,'YTick',1:90,'Fontsize',5)
set(gca,'YTickLabel',label90(end:-1:1,:))
ylim([0 91])
xlim([-1 1])
title('V_1 of BOLD Phase Coherence','Fontsize',12)
grid on

% subplot(4,4,[13 14 15])
% colormap(jet)
% FC_V=V_best'*V_best;  
% li=max(abs(FC_V(:)));
% imagesc(FC_V,[-li li])
% axis square
% title('FC pattern') 
% ylabel('Brain area #')
% xlabel('Brain area #')  

% 
% %% LIFETIMES
% 
% sigC=zeros(length(rangeK));
% 
% for k=1:length(rangeK)   
%     [Min_p_value(k), sigC(k)]=min(LT_pval(k,LT_pval(k,:)>0));    
% end
% 
% figure
% semilogy(rangeK,0.05*ones(1,length(rangeK)),'r--')
% hold on
% semilogy(rangeK,0.05./rangeK,'g--')
% semilogy(rangeK(1)-1:rangeK(end)+1,0.05./sum(rangeK)*ones(1,length(rangeK)+2),'b--','LineWidth',1)
% 
% 
% for k=1:length(rangeK) 
%     for c=1:rangeK(k)        
%         if LT_pval(k,c)>0.05
%             semilogy(rangeK(k),LT_pval(k,c),'*k');
%         end
%         if LT_pval(k,c)<0.05 && LT_pval(k,c)>(0.05/rangeK(k))
%             semilogy(rangeK(k),LT_pval(k,c),'*r');
%         end
%         if LT_pval(k,c)<(0.05/rangeK(k)) && LT_pval(k,c)>(0.05/sum(rangeK))
%             semilogy(rangeK(k),LT_pval(k,c),'*g');
%         end
%         if LT_pval(k,c)<=(0.05/sum(rangeK))
%             semilogy(rangeK(k),LT_pval(k,c),'*b');
%         end
%     end
% end
% 
% 
% ylabel('LT Patients vs Controls (p-value)')
% xlabel('Number of clusters K')
% box off
% 
% figure
% V_all=zeros(N_areas,length(rangeK));
% 
% for k=1:length(rangeK)
%     c=sigC(k);
%     V=Kmeans_results{k}.C(c,:);
%     V_all(:,k)=V;
%     
%     V=V(Order);
%     subplot(1,length(rangeK),k)
%     hold on
%     barh((V.*(V<=0)),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
%     barh((V.*(V>0)),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
%     ylim([0 N_areas+1])
%     xlim([-.15 .15])
%     
%     grid on
%     set(gca,'YTick',1:N_areas,'Fontsize',8)    
%     if k==1
%         set(gca,'YTickLabel',label90(end:-1:1,:))
%     else
%         set(gca,'YTickLabel',[])
%     end
%     if Min_p_value(k)<0.05 && Min_p_value(k)>(0.05/rangeK(k))
%         title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k))]},'Fontsize',10,'color','r')
%     elseif Min_p_value(k)>0.05
%         title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k))]},'Fontsize',10,'color','k')
%     elseif Min_p_value(k)<(0.05/k) && Min_p_value(k)>(0.05/sum(rangeK))
%         title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k))]},'Fontsize',10,'color','g')
%     elseif Min_p_value(k)<(0.05/sum(rangeK))
%         title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k))]},'Fontsize',10,'color','b')        
%     end
%     title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k))]},'Fontsize',10)
% end

% V_best=mean(V_all(:,Min_p_value<0.05),2)';
% 
% figure
% subplot(3,2,1)
% plot_nodes_in_cortex(V_best)
% rotate3d
% view(0,0)
% title('Mean signif net Lifetime')
% 
% subplot(3,2,3)
% plot_nodes_in_cortex(V_best)
% view(-90,0)
% rotate3d
% 
% subplot(3,2,[2 4 6])
% hold on
% barh(V.*(V<=0),'FaceColor','blue','EdgeColor','none')
% barh(V.*(V>0),'FaceColor','red','EdgeColor','none','Barwidth',.4)
% set(gca,'YTick',1:116,'Fontsize',7)
% set(gca,'YTickLabel',label90(end:-1:1,:))
% ylim([0 117])
% xlim([-.15 .15])
% title('relative BOLD phase','Fontsize',12)
% grid on
% 
% subplot(3,2,5)
% colormap(jet)
% FC_V=V_best'*V_best;  
% li=max(abs(FC_V(:)));
% imagesc(FC_V,[-li li])
% axis square
% title('FC pattern') 
% ylabel('Brain area #')
% xlabel('Brain area #')  
end