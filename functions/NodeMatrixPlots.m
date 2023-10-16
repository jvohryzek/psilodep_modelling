function [plotType] = NodeMatrixPlots( V )
%UNTITLED7 Summary of this function goes here

% Add necessary paths
addpath(genpath('/Volumes/Atilla/localadmin/matlab/Collaboration_Atasoy/LEiDA/plot_kit/code'));

% Set paths for surface and label files
path_surf_rh = '/Volumes/Atilla/localadmin/matlab/Collaboration_Atasoy/LEiDA/plot_kit/needed_data/bert_FREESURFER/surf/rh.pial';
path_surf_lh = '/Volumes/Atilla/localadmin/matlab/Collaboration_Atasoy/LEiDA/plot_kit/needed_data/bert_FREESURFER/surf/lh.pial';
% Set here which parcellation!
% scale 1 -> _36.annot
% scale 2 -> _60.annot
% scale 3 -> _125.annot
% scale 4 -Y _250.annot
scale = 1;
path_annot_rh = '/Volumes/Atilla/localadmin/matlab/Collaboration_Atasoy/LEiDA/plot_kit/needed_data/bert_FREESURFER/label/rh.myaparc_36.annot';
path_annot_lh = '/Volumes/Atilla/localadmin/matlab/Collaboration_Atasoy/LEiDA/plot_kit/needed_data/bert_FREESURFER/label/lh.myaparc_36.annot';

% Load necessary data
load('/Volumes/Atilla/localadmin/matlab/Collaboration_Atasoy/LEiDA/plot_kit/needed_data/labels_index_CORTICAL_Laus2008_all_scales.mat');


%% plot num. 1
map = parula;
nodeAff = V(1,:);
nodeAff = nodeAff(ixc{scale});

CM = zeros(length(ixc{scale}),3); % NOTE that only cortical nodes will be plotted. It is not possible to visualize subcortical nuclei.
CM(1:68,:) = mapsc2rgb(rand(size(ixc{scale})), map); % This function map a given function test_function to a selected colormap map
%
plotType = 3;
matrixPreTh = zeros(68,68); % just a dummy for the edges
colorsurf_2hemi_node(path_surf_rh, path_annot_rh, path_surf_lh, path_annot_lh,CM, llist{scale}, matrixPreTh, plotType, nodeAff);


% xlabel('ROIs','FontSize',26);ylabel('Degree','FontSize',26);title(strcat('Degree distribution of ',modality),'FontSize',26)

end

