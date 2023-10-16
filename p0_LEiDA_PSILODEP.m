
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEADING EIGENVECTOR DYNAMICS ANALYSIS (LEiDA)
%
%
% Jakub Vohryzek jakub.vohryzek@upf.edu
% original: Joana Cabral May 2016 
% joana.cabral@psych.ox.ac.uk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% Part 0. INITIALISATION %%%%%%%%%%%%%%%%%%%%%%%%%%
%% Decision Tree
figure_save             = 0;  % 0 if no saving figures and 1 yes saving figures
result_save             = 0;  % 0 if no saving results and 1 yes saving results
result2model_save       = 0;  % 0 if no saving results and 1 yes saving results 
compute_kmeans          = 0;  % 0 if no saving load results and 1 yes run kmeans
compute_stats           = 0;  % 0 if no saving load results and 1 yes run stats
%% load path
Directory               ='/Users/jakub/Matlab/Collaboration_Cabral';
inputData               = '/Users/jakub/Datasets';
namesRest_before        = dir([inputData '/PSILODEP/MAT/S*rest*before*']);
namesRest_after         = dir([inputData '/PSILODEP/MAT/S*rest*after*']);
namesMusic_before       = dir([inputData '/PSILODEP/MAT/S*music*before*']);
namesMusic_after        = dir([inputData '/PSILODEP/MAT/S*music*after*']);
SaveFile                = [Directory, '/LEiDA/LEiDA-Vohryzek_PSILODEP/Results/24_20_2019_detrend_thesis/'];

namesRest               = [namesRest_before; namesRest_after];
namesMusic              = [namesMusic_before; namesMusic_after];
names                   = [namesRest; namesMusic];
%% AAL labels
load(strcat(Directory,'/LEiDA/LEiDA-Vohryzek_PSILODEP/Data/AAL_labels.mat'))

%% FILTER SETTINGS 
TR                      = 2;                 
fnq                     = 1/(2 * TR);               % Nyquist frequency
flp                     = 0.04; %                    % lowpass frequency of filter
fhi                     = 0.07;%                    % highpass
Wn                      = [flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k                       = 2;                          % 2nd order butterworth filter
[bfilt,afilt]           = butter(k,Wn);   % construct the filter

%% Initialising

numAreas                = 90;                  % number of brain regions
numTp                   = 231;                    % number of timepoints = 237 points
excTp                   = 3; % number of timepoints to exclude from before and after
numSbj                  = length(namesRest_before);% number of subjects
colorVal                = {'r';'m';'y';'g';'b';'c';'r';'m';'y';'g';'b';'c';'r';'m';'y';'g';'b';'c';'r';'m';'y';'g';'b';'c'};
font_size               = 24;
%% REMISSION
% version: in 5 week (QIDS)
responders             = [1 1 0 0 0 1 0 1 1 1 0 0 0 0 0];

respNum                 = sum(responders);
norespNum = sum(responders == 0);
%% evaluating QIDS
QIDS_baseline = [19, 28, 18, 18, 25, 22, 17, 26, 28, 26, 29, 36,28, 24, 28];
QIDS_baseline_rp_mean = mean(QIDS_baseline(responders ==1));
QIDS_baseline_rp_std = std(QIDS_baseline(responders ==1));

QIDS_baseline_nrp_mean = mean(QIDS_baseline(responders ==0));
QIDS_baseline_nrp_std = std(QIDS_baseline(responders ==0));

QID_baseline_diff = permutation_htesting_np([QIDS_baseline(responders ==1),...
                    QIDS_baseline(responders == 0)],...
                    [ones(1,numel(QIDS_baseline(responders ==1))) 2*ones(1,numel(QIDS_baseline(responders ==0)))],5000,0.05,'ranksum');
QID_baseline_diff_pval = min(QID_baseline_diff.pvals);
%% %%%%%%%%%%%%%%%%% kmeans to be adjusted accordingly %%%%%%%%%%%%%%%%%%%%
maxK                    = 10;
replicates              = 20;
distance                = 'sqeuclidean';
Condition               = {'Rest'}; 
workDate                = '24_20_2019_detrend_thesis'

% save fileCluster
SaveFileCluster         = strcat('/Users/jakub/Matlab/Collaboration_Cabral/LEiDA/LEiDA-Vohryzek_PSILODEP/Results/',workDate,'/Clusters/');
if ~exist(SaveFileCluster,'dir')
    mkdir(SaveFileCluster)
    mkdir(strcat(SaveFileCluster,distance))
end

%% 

BOLD_processed_idx      = [] ; Phase_BOLD_idx    = []; staticFC_idx = [];
Leading_Eig_idx         = [] ; Var_Eig_idx       = []; iFC_idx = [];
OP_idx                  = [] ; Synchro_idx       = []; Metasta_idx = [];
% To reorder matrix plots
order                   = [1:2:numAreas numAreas:-2:2]; % AAL is interleaving right and left brain regions
% conditions
numCond                 = length(Condition);
if length(Condition) == 2
    numTask             = 4;                            % Restbefore, Restafter, Musicbefore, Musicafter
else
    numTask             = 2;                            % Restbefore, Restafter
end
numIdx                  = numTask*numSbj;               % number of SubjectsxTask

%% %%%%%%%%%%%%%%%%%%% Part 1. DATA %%%%%%%%%%%%%%%%%%%%%%%%%%
data = struct;
data.remission                      = struct;      data.no_remission                  = struct;
data.remission.beforeRest           = struct;      data.remission.afterRest           = struct;
data.no_remission.beforeRest        = struct;      data.no_remission.afterRest        = struct;

data.remission.beforeRest.LEiall    = [];          data.remission.afterRest.LEiall    = [];
data.no_remission.beforeRest.LEiall = [];          data.no_remission.afterRest.LEiall = [];

XBOLD                               = cell(numTask,numSbj); % for the optimisation part
for task = 1:numTask          % n_Task: Restbefore,Restafter

    sbjrem = 0; sbjnorem = 0;
    for sbj = 1:numSbj % n_Subjects
        % Indexing
        idx = numSbj*(task-1)+sbj; % indexing to save
        idx

        % load
        BOLD                    = load([names(idx).folder '/' names(idx).name])';
        XBOLD{task,sbj}         = BOLD;
        
        % BOLD signal to BOLD phase using the Hilbert transform
        [~, Phase_BOLD, ~]      = BOLD2hilbert(BOLD.ts, numAreas, bfilt, afilt,excTp);
        [~, ~, Leading_Eig ]    = instFC( Phase_BOLD, numAreas, numTp,'default','halfSwitching');
        
        %% SAVING INTO DATA STRUCT
        if responders(sbj)
                sbjrem = sbjrem+1;
            if task == 1
                data.remission.beforeRest.LEiall    = cat(2,data.remission.beforeRest.LEiall,Leading_Eig);
            else
                data.remission.afterRest.LEiall     = cat(2,data.remission.afterRest.LEiall,Leading_Eig);
            end
        else
            sbjnorem = sbjnorem+1;
            if task == 1
                data.no_remission.beforeRest.LEiall = cat(2,data.no_remission.beforeRest.LEiall,Leading_Eig);
            else
                data.no_remission.afterRest.LEiall  = cat(2,data.no_remission.afterRest.LEiall,Leading_Eig);
            end
        end
    end
end

%% %%%%%%%%%%%%%%%%%%% Part 2. CLUSTERS %%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%% KMEANS ALGORITHM %%%
% distance = 'sqeuclidean';
if compute_kmeans
    Leading_Eig_idx = [];
    % SBJ: remissionBF,remissionAF,noremissionBF, noremissionAF
    Leading_Eig_idx = [data.remission.beforeRest.LEiall data.remission.afterRest.LEiall data.no_remission.beforeRest.LEiall data.no_remission.afterRest.LEiall];
    % only if you want to run consistency across several runs    
    for it=1:10
        if strcmp(distance, 'sqeuclidean')
        [Kmeans_results, distM_FCD, ind_max, dunnScore(:,it), evaDB(it,:), evaSilh(it,:)] = kmeansClusters(Leading_Eig_idx, maxK, replicates, distance ); % ,'Display','final'); %,'Options',opt); 
        elseif strcmp(distance, 'cosine')
        [Kmeans_results, distM_FCD, ind_max, dunnScore(:,it), ~, evaSilh(it,:)] = kmeansClusters(Leading_Eig_idx, maxK, replicates, distance ); % ,'Display','final'); %,'Options',opt); 
        end
    end
    if result_save
        save(strcat(Directory, '/LEiDA/LEiDA-Vohryzek_PSILODEP/Results/'...
            ,'24_20_2019_detrend','/Clusters/',distance,...
            '/PSILODEPresults_AALClusters_',distance)) 
    end
else
    load(strcat(Directory, '/LEiDA/LEiDA-Vohryzek_PSILODEP/Results/'...
        ,'24_20_2019_detrend','/Clusters/',distance,...
        '/PSILODEPresults_AALClusters_',distance))  
end

%% %%% PLOTTING KMEANS EVALUATION RESULTS %%%
figure
subplot(131)
tmp_median = nanmedian(dunnScore(2:10,:)'); tmp_25qntile = quantile(dunnScore(2:10,:)',0.25); tmp_75qntile = quantile(dunnScore(2:10,:)',0.75);
xVec = 1:size(tmp_median,2); xConfresp= [xVec xVec(end:-1:1)]; yConfresp = [(tmp_25qntile) (tmp_75qntile(end:-1:1))];
p3=fill(xConfresp,yConfresp,'red'); p3.FaceColor = [0 0 0]; p3.EdgeColor = 'none'; p3.FaceAlpha = 0.5;
hold on
plot(tmp_median,'Color',[0 0 0],'LineWidth',2)
legend('IQL','median');title('Dunns Score'),xlabel('# of clusters'),ylabel('Dunns Score');
set(gca, 'XTick',1:9, 'XTickLabels',2:10,'FontSize',18);xlim([1,9]);
subplot(132)
tmp_median = nanmedian(evaDB);tmp_25qntile = quantile(evaDB,0.25);tmp_75qntile = quantile(evaDB,0.75);
xVec = 1:size(tmp_median,2);xConfresp= [xVec xVec(end:-1:1)];yConfresp = [(tmp_25qntile) (tmp_75qntile(end:-1:1))];
p3=fill(xConfresp,yConfresp,'red');p3.FaceColor = [0 0 0];p3.EdgeColor = 'none';p3.FaceAlpha = 0.5;
hold on
plot(tmp_median,'Color',[0 0 0],'LineWidth',2)
legend('IQL','median');title('Davies-Bouldin Score'),xlabel('# of clusters'),ylabel('Davies-Bouldin Score');
set(gca, 'XTick',1:9, 'XTickLabels',2:10,'FontSize',18);xlim([1,9]);
subplot(133)
tmp_median = nanmedian(evaSilh);tmp_25qntile = quantile(evaSilh,0.25);tmp_75qntile = quantile(evaSilh,0.75);
xVec = 1:size(tmp_median,2);xConfresp= [xVec xVec(end:-1:1)];yConfresp = [(tmp_25qntile) (tmp_75qntile(end:-1:1))];
p3=fill(xConfresp,yConfresp,'red');p3.FaceColor = [0 0 0];p3.EdgeColor = 'none';p3.FaceAlpha = 0.5;
hold on
plot(tmp_median,'Color',[0 0 0],'LineWidth',2)
legend('IQL','median');title('Silhouette Score'),xlabel('# of clusters'),ylabel('Silhouette Score');
set(gca, 'XTick',1:9, 'XTickLabels',2:10,'FontSize',18);xlim([1,9]);

%% %%%%%%%%%%%%%%%%%%%% Part 3. ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vsorted = {};
ind = {};numClust = {};
for m = 3:maxK
    % SORTING CLUSTERS
    Clusters{m} = Kmeans_results{m};
    typeSort = 'descend';
    [numClust{m}, numAreas] = size(Clusters{m}.C);
    [clustExpressionOrg, ind{m}, Vsorted{m}] = kmeansClustersSorted(Clusters{m},numClust{m},typeSort);
end
%% plot the brain states renderings
PSILODEP_Plot_repertoire_RL_LR_Pub(Vsorted{3})

%%
if compute_stats

    for m = 2:10

        % 3. FRACTIONAL OCCUPANCY - permutation testing
        [Cluster_Param{m}] = clusterStatsFA(Clusters{m},ind{m},numSbj,numTask,numTp,numClust{m},respNum,norespNum,'cross-sectional');
        for i = 1:m
            % non-parametric
            disp('Compare Cluster Probability of before All with after All')
            disp(['Now running for ' num2str(i) ' cluster'])
            NPstatsRest = permutation_htesting_np([Cluster_Param{m}{1}.foRestBefore(i,:),Cluster_Param{m}{1}.foRestAfter(i,:)],...
                [ones(1,numel(Cluster_Param{m}{1}.foRestBefore(i,:))) 2*ones(1,numel(Cluster_Param{m}{1}.foRestAfter(i,:)))],1000,0.05,'signrank');
            Cluster_Statsv2.pNPRest(m,i) = min(NPstatsRest.pvals);

            disp('Compare Cluster Probability of before Resp with after Resp')
            disp(['Now running for ' num2str(i) ' cluster'])
            NPstatsRestResp = permutation_htesting_np([Cluster_Param{m}{1}.foRestAfterResp(i,:),Cluster_Param{m}{1}.foRestBeforeResp(i,:)],...
                [ones(1,numel(Cluster_Param{m}{1}.foRestAfterResp(i,:))) 2*ones(1,numel(Cluster_Param{m}{1}.foRestBeforeResp(i,:)))],1000,0.05,'signrank');
            Cluster_Statsv2.pNPRestResp(m,i) = min(NPstatsRestResp.pvals);

            disp('Compare Cluster Probability of before NoResp with after NoResp')
            disp(['Now running for ' num2str(i) ' cluster'])
            NPstatsRestNoResp = permutation_htesting_np([Cluster_Param{m}{1}.foRestAfterNoResp(i,:),Cluster_Param{m}{1}.foRestBeforeNoResp(i,:)],...
                [ones(1,numel(Cluster_Param{m}{1}.foRestAfterNoResp(i,:))) 2*ones(1,numel(Cluster_Param{m}{1}.foRestBeforeNoResp(i,:)))],1000,0.05,'signrank');
            Cluster_Statsv2.pNPRestNoResp(m,i) = min(NPstatsRestNoResp.pvals);

            disp('Compare Cluster Probability of Resp Diff with NoResp Diff')
            disp(['Now running for ' num2str(i) ' cluster'])
            NPstatsRestdiff = permutation_htesting_np([Cluster_Param{m}{1}.foRestDiffResp(i,:),Cluster_Param{m}{1}.foRestDiffNoResp(i,:)],...
                [ones(1,numel(Cluster_Param{m}{1}.foRestDiffResp(i,:))) 2*ones(1,numel(Cluster_Param{m}{1}.foRestDiffNoResp(i,:)))],5000,0.05,'ranksum');
            Cluster_Statsv2.pNPRestdiff(m,i) = min(NPstatsRestdiff.pvals);

            disp('Compare Cluster Probability of after Resp with after NoResp')
            disp(['Now running for ' num2str(i) ' cluster'])
            NPstatsRestAfRestNoRest = permutation_htesting_np([Cluster_Param{m}{1}.foRestAfterResp(i,:),Cluster_Param{m}{1}.foRestAfterNoResp(i,:)],...
                [ones(1,numel(Cluster_Param{m}{1}.foRestAfterResp(i,:))) 2*ones(1,numel(Cluster_Param{m}{1}.foRestAfterNoResp(i,:)))],1000,0.05,'ranksum');
            Cluster_Statsv2.pNPRestAfRestNoRest(m,i) = min(NPstatsRestAfRestNoRest.pvals);

            disp('Compare Cluster Probability of After NoResp with before All')
            disp(['Now running for ' num2str(i) ' cluster'])
            NPstatsRestBfAllAfRest = permutation_htesting_np([Cluster_Param{m}{1}.foRestAfterResp(i,:),Cluster_Param{m}{1}.foRestBefore(i,:)],...
                [ones(1,numel(Cluster_Param{m}{1}.foRestAfterResp(i,:))) 2*ones(1,numel(Cluster_Param{m}{1}.foRestBefore(i,:)))],1000,0.05,'ranksum');
            Cluster_Statsv2.pNPRestBfAllAfRest(m,i) = min(NPstatsRestBfAllAfRest.pvals);

            disp('Compare Cluster Probability of before Resp with before NoResp')
            disp(['Now running for ' num2str(i) ' cluster'])
            NPstatsRestBfRestNoRest = permutation_htesting_np([Cluster_Param{m}{1}.foRestBeforeResp(i,:),Cluster_Param{m}{1}.foRestBeforeNoResp(i,:)],...
                [ones(1,numel(Cluster_Param{m}{1}.foRestBeforeResp(i,:))) 2*ones(1,numel(Cluster_Param{m}{1}.foRestBeforeNoResp(i,:)))],1000,0.05,'ranksum');
            Cluster_Statsv2.pNPRestBfRestNoRest(m,i) = min(NPstatsRestBfRestNoRest.pvals);

        end  
    end
    if result_save
        save(strcat(SaveFile,workDate,'/Clusters/',distance,'/','/PSILODEPresults_AALClusters_pvaluesv2_',distance), 'Cluster_Statsv2', 'Cluster_Param','Vsorted','Cluster_Statsv3')
    end
else
    load(strcat(Directory, '/LEiDA/LEiDA-Vohryzek_PSILODEP/Results/'...
    ,'24_20_2019_detrend','/Clusters/',distance,...
    '/PSILODEPresults_AALClusters_pvaluesv2_',distance)) 
end

%% h1pNP = figure
load('/Users/jakub/Matlab/Collaboration_Cabral/LEiDA/LEiDA-Vohryzek_PSILODEP/Results/24_20_2019_detrend/Clusters/sqeuclidean/PSILODEPresults_AALClusters_pvaluesv2_sqeuclidean.mat')
%% non-parametric statistics
h1np = figure

subplot(2,4,1); Plot_p_values_FunHarD(Cluster_Statsv2.pNPRest,2:10,'Probability');title('Probability of Patients Before and After');
subplot(2,4,2); Plot_p_values_FunHarD(Cluster_Statsv2.pNPRestResp,2:10,'Probability');title('Probability of Responders Before and After');
subplot(2,4,3); Plot_p_values_FunHarD(Cluster_Statsv2.pNPRestNoResp,2:10,'Probability');title('Probability of Non-Responders Before and After');
subplot(2,4,4); Plot_p_values_FunHarD(Cluster_Statsv2.pNPRestdiff,2:10,'Probability');title('Probability of Resp and NoResp Diff');
subplot(2,4,5); Plot_p_values_FunHarD(Cluster_Statsv2.pNPRestAfRestNoRest,2:10,'Probability');title('Probability of after Resp and NoResp');
subplot(2,4,6); Plot_p_values_FunHarD(Cluster_Statsv2.pNPRestBfRestNoRest,2:10,'Probability');title('Probability of before Resp and NoResp');
subplot(2,4,7); Plot_p_values_FunHarD(Cluster_Statsv2.pNPRestBfAllAfRest,2:10,'Probability');title('Probability of After Resp with before All');
%% 
%% this part is hard-coded for solution of k=3
mainNames = {'Before Responders','After Responders','Before Non-Responders','After Non-Responders'}
cMap =  [0 0.2 0; 0 0.2 0; 0.4 0 0.2; 0.4 0 0.2];

figure
for i = 1:3
    Cluster_Param_plot.BeforeResp   = Cluster_Param{3}{1}.foRestBeforeResp(i,:);
    Cluster_Param_plot.AfterResp    = Cluster_Param{3}{1}.foRestAfterResp(i,:);
    Cluster_Param_plot.BeforeNoResp = Cluster_Param{3}{1}.foRestBeforeNoResp(i,:);
    Cluster_Param_plot.AfterNoResp  = Cluster_Param{3}{1}.foRestAfterNoResp(i,:);
    Cluster_Param_plot.BeforeResp   = Cluster_Param{3}{1}.foRestBeforeResp(i,:);
    Cluster_Param_stat(1,:)         = Cluster_Statsv2.pNPRestResp(3,1:3);
    Cluster_Param_stat(2,:)         = Cluster_Statsv2.pNPRestNoResp(3,1:3);
    Cluster_Param_stat(3,:)         = Cluster_Statsv2.pNPRestAfRestNoRest(3,1:3);

    subplot(1,3,i)
    PMS_states_plot(Cluster_Param_plot,Cluster_Param_stat(:,i),mainNames,cMap);ylim([0 1])
end
% additional paired-boxplot
% [h4pairedboxplot] = plotClustersFObars(Cluster_Param,3,respNum,norespNum,numClust{3},colorVal,'pairedboxplot');
%% KOP stats based on 'data' structure
load('/Users/jakub/Matlab/Collaboration_Cabral/LEiDA/LEiDA-Vohryzek_PSILODEP/Results/24_20_2019_detrend/Clusters/sqeuclidean/PSILODEPresults_AALDatarem.mat', 'data')
%%
figure
subplot(232)
KOP_Synchrony_plot.BeforeResp   = data.remission.beforeRest.Synchrony;
KOP_Synchrony_plot.AfterResp    = data.remission.afterRest.Synchrony;
KOP_Synchrony_plot.BeforeNoResp = data.no_remission.beforeRest.Synchrony;
KOP_Synchrony_plot.AfterNoResp  = data.no_remission.afterRest.Synchrony;
KOP_stat(1,:)         = Cluster_Statsv2.pNPRestResp(3,1:3);
KOP_stat(2,:)         = Cluster_Statsv2.pNPRestNoResp(3,1:3);
KOP_stat(3,:)         = Cluster_Statsv2.pNPRestAfRestNoRest(3,1:3);

KOP_states_plot(KOP_Synchrony_plot,KOP_Synchrony_plot,mainNames,cMap);ylim([0 .7])
ylabel('Synchrony')
subplot(233)
KOP_Metastability_plot.BeforeResp   = data.remission.beforeRest.Metastability;
KOP_Metastability_plot.AfterResp    = data.remission.afterRest.Metastability;
KOP_Metastability_plot.BeforeNoResp = data.no_remission.beforeRest.Metastability;
KOP_Metastability_plot.AfterNoResp  = data.no_remission.afterRest.Metastability;
KOP_stat(1,:)         = Cluster_Statsv2.pNPRestResp(3,1:3);
KOP_stat(2,:)         = Cluster_Statsv2.pNPRestNoResp(3,1:3);
KOP_stat(3,:)         = Cluster_Statsv2.pNPRestAfRestNoRest(3,1:3);

KOP_states_plot(KOP_Metastability_plot,KOP_Synchrony_plot,mainNames,cMap);ylim([0 .3])
ylabel('Metastability')

subplot(2,3,[4,5,6])
histogram(data.remission.beforeRest.histFCD,1000/8,'Normalization','probability','FaceColor',cMap(1,:),'EdgeColor','none')
hold on
histogram(data.remission.afterRest.histFCD,1000/8,'Normalization','probability','FaceColor',cMap(2,:),'EdgeColor','none')
hold on
histogram(data.no_remission.beforeRest.histFCD,1000/8,'Normalization','probability','FaceColor',cMap(3,:),'EdgeColor','none')
hold on
histogram(data.no_remission.afterRest.histFCD,1000/8,'Normalization','probability','FaceColor',cMap(4,:),'EdgeColor','none')
set(gca,'FontSize',24);xlim([-1 1])
legend(mainNames,'Location','southwest');
%set(hAdd0,'Units', 'Inches', 'Position', [0, 0, 55, 55], 'PaperUnits', 'Inches', 'PaperSize', [55, 55])
xlabel('FCD values');ylabel('Probability');
%% FCD statistics - comparing the mean of each subjects distribution
for i=1:6
   FCD_resp_before_mean(i) = mean(nonzeros(triu(data.remission.beforeRest.FCD(:,:,i))));
   FCD_resp_after_mean(i)  = mean(nonzeros(triu(data.remission.afterRest.FCD(:,:,i))));
   FCD_resp_before_std(i) = std(nonzeros(triu(data.remission.beforeRest.FCD(:,:,i))),[],1);
   FCD_resp_after_std(i)  = std(nonzeros(triu(data.remission.afterRest.FCD(:,:,i))),[],1);

end
for i=1:9
   FCD_noresp_before_mean(i) = mean(nonzeros(triu(data.no_remission.beforeRest.FCD(:,:,i))));
   FCD_noresp_after_mean(i)  = mean(nonzeros(triu(data.no_remission.afterRest.FCD(:,:,i))));
   FCD_noresp_before_std(i) = std(nonzeros(triu(data.no_remission.beforeRest.FCD(:,:,i))),[],1);
   FCD_noresp_after_std(i)  = std(nonzeros(triu(data.no_remission.afterRest.FCD(:,:,i))),[],1);

end
% non-parametric
NPstats_FCD_resp_before_after = permutation_htesting_np([FCD_resp_before_mean,FCD_resp_after_mean],...
    [ones(1,numel(FCD_resp_before_mean)) 2*ones(1,numel(FCD_resp_after_mean))],1000,0.05,'signrank');
pNPstats_FCD_resp_before_after = min(NPstats_FCD_resp_before_after.pvals);

NPstats_FCD_resp_before_nonresp_before = permutation_htesting_np([FCD_resp_before_mean,FCD_noresp_before_mean],...
    [ones(1,numel(FCD_resp_before_mean)) 2*ones(1,numel(FCD_noresp_before_mean))],1000,0.05,'ranksum');
pNPstats_FCD_resp_before_nonresp_before = min(NPstats_FCD_resp_before_nonresp_before.pvals);

NPstats_FCD_resp_before_nonresp_after = permutation_htesting_np([FCD_resp_before_mean,FCD_noresp_after_mean],...
    [ones(1,numel(FCD_resp_before_mean)) 2*ones(1,numel(FCD_noresp_after_mean))],1000,0.05,'ranksum');
pNPstats_FCD_resp_before_nonresp_after = min(NPstats_FCD_resp_before_nonresp_after.pvals);

NPstats_FCD_resp_after_nonresp_after = permutation_htesting_np([FCD_resp_after_mean,FCD_noresp_after_mean],...
    [ones(1,numel(FCD_resp_after_mean)) 2*ones(1,numel(FCD_noresp_after_mean))],1000,0.05,'ranksum');
pNPstats_FCD_resp_after_nonresp_after = min(NPstats_FCD_resp_after_nonresp_after.pvals);

NPstats_FCD_resp_after_nonresp_before= permutation_htesting_np([FCD_resp_after_mean,FCD_noresp_before_mean],...
    [ones(1,numel(FCD_resp_after_mean)) 2*ones(1,numel(FCD_noresp_before_mean))],1000,0.05,'ranksum');
pNPstats_FCD_resp_after_nonresp_before = min(NPstats_FCD_resp_after_nonresp_before.pvals);


%%
% Global Brain Connectivity
for i=1:6
    GBC.BeforeResp(i)   = mean(sum(mean(data.remission.beforeRest.iFC{i},3)));
    GBC.AfterResp(i) = mean(sum(mean(data.remission.afterRest.iFC{i},3)));
end
for i=1:9
    GBC.BeforeNoResp(i)   = mean(sum(mean(data.no_remission.beforeRest.iFC{i},3)));
    GBC.AfterNoResp(i) = mean(sum(mean(data.no_remission.afterRest.iFC{i},3)));
end
subplot(231)

GBC_stat(1,:)         = Cluster_Statsv2.pNPRestResp(3,1:3);
GBC_stat(2,:)         = Cluster_Statsv2.pNPRestNoResp(3,1:3);
GBC_stat(3,:)         = Cluster_Statsv2.pNPRestAfRestNoRest(3,1:3);
KOP_states_plot(GBC,KOP_Synchrony_plot,mainNames,cMap);%ylim([0 .7])
ylabel('GBC')
%% SAVING FOR THE MODELLING PART
if result2model_save
    % for cluster solution 3

    P1emp = Cluster_Param{1,3}{1,1}.foRestBefore;
    P2emp = Cluster_Param{1,3}{1,1}.foRestAfter;
    PTR1emp = Cluster_Param{1,3}{1,1}.tmRestBefore;
    PTR2emp = Cluster_Param{1,3}{1,1}.tmRestAfter;
    Vemp = Vsorted{1,3};
    save(strcat(SaveFileCluster,distance,'/empiricalLEiDA_PSILODEP_filt0_04_0_07v3.mat'),...
        'XBOLD', 'Vemp', 'P1emp', 'P2emp', 'PTR1emp', 'PTR2emp', 'numClust');
end