%% %%%%%%%%%%%%%%%%%
% A script to compute overlap of the PSILODEP perturbation results with
%
% Jakub Vohryzek jakub.vohryzek@upf.edu
%
%% %%%%%%%%%%%%%%%%

%% Define the Yeo Networks in the new Parcellation
Directory = '/Users/jakub/Matlab/Collaboration_Kringelbach/Psilodep_modelling/';
%%  load the group differences 
load([Directory, 'perturbation_difference_results_31_03_2021.mat'])
x = abs(a_perturb_stats.diff_symm');

%%  load the 5HT2A maps
load([Directory 'mean5T_all.mat'])

%% 5HT2A map correlation

figure
subplot(1,2,1)

y1 = symm_mean5HT2A; 
% correlation
[RHO,PVAL] = corr(x,y1,'Type','Spearman');

p = polyfit(x,y1,1);
% Evaluate the fitted polynomial p and plot:
f = polyval(p,x);
plot(x,y1,'.',x,f,'-','MarkerSize',20,'LineWidth',4);axis square
xlabel({'(KLdiv Resp/ - KLdiv Non-Resp.)','at optimal a'})
ylabel('5HT2a density')
ylabel('5HT2a')

title({['Spearman Correlation of ' num2str(round(RHO,3))],['p-value of ' num2str(round(PVAL,3))]})
set(gca, 'FontSize',14);xlim([min(x) max(x)]);
%% 5HT1A map correlation
subplot(1,2,2)
%  load the 5HT1A maps
y2 = symm_mean5HT1A;
% correlation
[RHO,PVAL] = corr(x,y2,'Type','Spearman');

p = polyfit(x,y2,1);
% Evaluate the fitted polynomial p and plot:
f = polyval(p,x);
plot(x,y2,'.',x,f,'-','MarkerSize',20,'LineWidth',4);axis square
xlabel({'(KLdiv Resp/ - KLdiv Non-Resp.)','at optimal a'});
ylabel('5HT1a density')
title('5HT1a')
title({['Spearman Correlation of ' num2str(round(RHO,3))],['p-value of ' num2str(round(PVAL,3))]})
set(gca, 'FontSize',14);xlim([min(x) max(x)]);
%% plotting the perturbation difference on th brain
figure
rendersurface_aal90_gifti(x,1,min(x),max(x),0,'YlOrRd9',2)
colorbar;
%% plotting the first 10 ROIs
[~, I] = sort(x,'descend'); %
x_10_highest = (x >= x(I(20))); % I(20)cut-off value for the first 10 ROIs
figure
rendersurface_aal90_gifti(x_10_highest,1,min(x),max(x),0,'YlOrRd9',2)
%%%%%%%%%%%%%%%%%%%%%%%%%% RSNs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the Yeo Networks in the new Parcellation
N_areas                          = 90;
% First load the mask of the chosen parcellation in MNI space 2mm
Parcellation                     = 'AAL116';
V_Parcels                        = struct2array(load('/Users/jakub/Matlab/Collaboration_Cabral/LEiDA/Modelling-Vohryzek_PSILODEP/Vohryzek_Code/psilodep_braincomm_submission/ParcelsMNI2mm.mat',['V_' Parcellation]));
V_Parcels((V_Parcels > N_areas)) = 0;
%% load the mask of the Yeo parcellation in MNI space 2mm
V_Yeo                            = struct2array(load('ParcelsMNI2mm','V_Yeo7'));
%% prepare the inputs
order                            = [1:2:90 90:-2:1]; unique(order);
[~, idx_symm]                      = unique(order);
%%

Centroids{1}.C = x(idx_symm)'; %

Centroids{1}.C.*(Centroids{1}.C > 0.0854); % first 5 homological networks

%% additional receptors results

figure
subplot(1,3,1)
%  load the 5HT2A maps
y3 = symm_mean5HT1B;
% correlation
[RHO,PVAL] = corr(x,y3,'Type','Spearman');

p = polyfit(x,y3,1);
f = polyval(p,x);
plot(x,y3,'.',x,f,'-','MarkerSize',20,'LineWidth',4);axis square
xlabel({'(KLdiv Resp/ - KLdiv Non-Resp.)','at optimal a'})
ylabel('5HT1b density')
title('5HT1b')
set(gca, 'FontSize',14)
%%

subplot(1,3,2)
%  load the 5HT2A maps
y4 = symm_mean5HT4; 
% correlation
[RHO,PVAL] = corr(x,y4,'Type','Spearman');

p = polyfit(x,y4,1);
f = polyval(p,x);
plot(x,y4,'.',x,f,'-','MarkerSize',20,'LineWidth',4);axis square
xlabel({'(KLdiv Resp/ - KLdiv Non-Resp.)','at optimal a'})
ylabel('5HT4 density')
title('5HT4')
set(gca, 'FontSize',14)
%%
subplot(1,3,3)
%  load the 5HT2A maps
y5 = symm_mean5HTT; 
% correlation
[RHO,PVAL] = corr(x,y5,'Type','Spearman');

p = polyfit(x,y5,1);
f = polyval(p,x);
plot(x,y5,'.',x,f,'-','MarkerSize',20,'LineWidth',4);axis square
xlabel({'(KLdiv Resp/ - KLdiv Non-Resp.)','at optimal a'})
ylabel('5HTT density')
title('5HTT')
set(gca, 'FontSize',14)