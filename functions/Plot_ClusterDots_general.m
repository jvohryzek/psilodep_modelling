function Plot_ClusterDots_general

% Code to visualize all LEiDA observations colored according to the k-means clustering solution
% - plots each observation as a dot in 2 Principal Dimensions
% - plots the cluster centroids as a +
% - plots the corresponding matrices below
% - plots the probability of each
% 
% Joana Cabral and L.D. Lord, October 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose k
k = 5;

% Load the Cluster Centroids and sort them
%load LEiDA_psilo_newkresults.mat Kmeans_results
%load LEiDA_psilo_data.mat Time_sessions
IDX_Psi= Kmeans_results{k}.IDX% Kmeans_results{k}.IDX;

% Sort states according to their probability of occurrence in condition 1
condition=1;
for c=1:k
    ProbC(c)=mean(IDX_Psi(Time_sessions(1,:)==condition)==c);
end
[~, ind_sort]=sort(ProbC,'descend'); 
Centers= Kmeans_results{k}.C(ind_sort,:) % Kmeans_results{k}.C(ind_sort,:);
clear Kmeans_results

Order=[1:2:89 90:-2:2];

% DEFINE HERE THE COLORS OF THE CLUSTERS IN RGB CODE
cmap=[0 0 1; 0 1 1 ; 0 1 0 ; 1 1 0;  1 0 0; 1 0 1; .7 .7 .7 ; 1 0.5 0 ];


%% FIGURE OF Cluster Clouds

% NOW LOAD ALL THE EIGENVECTORS FROM ALL SESSIONS
load LEiDA_psilo_data.mat Leading_Eig Time_sessions
% Time sessions contains the time of each condition

% CONCATENATE ALL VECTORS WITH CENTROIDS AT THE END
V1_psi=[Leading_Eig; Centers];

% OBTAIN THE FIRST 3 PC OF THE COVARIANCE MATRIX 
Var=cov(V1_psi');
[pc3, ~]=eigs(Var,3);
% pc3 is a vector of 3x(TOTAL_TIME+K), where each row represent the
% coordinates X, Y and Z of each observation. 
% The last K values in pc3 are the centroids

% PLOT DOTS IN 3D IN BLACK Before Clustering
subplot(1,2,1)
hold on
for t=find(Time_sessions(1,:) == 1)  %plot for first condition
    plot3(pc3(t,3),pc3(t,2),pc3(t,1),'.','Markersize',10,'Color','k')
end
zlabel('1st PC')
ylabel('2nd PC') 
xlabel('3nd PC')
zlim([ -.06 0.06])
%set(gca,'Color','k')
title('Eigen Space')
view(85,10)

% PLOT DOTS IN 3D IN COLOR After Clustering
subplot(1,2,2)
hold on

% PLOT DOTS IN 3D WITH CORRESPONDING COLOR
for t=find(Time_sessions(1,:) == 1)  %plot for first condition
    plot3(pc3(t,3),pc3(t,2),pc3(t,1),'.','Markersize',10,'Color',cmap(ind_sort==IDX_Psi(t),:))
end
zlabel('1st PC')
ylabel('2nd PC') 
xlabel('3nd PC')
zlim([ -.06 0.06])
%set(gca,'Color','k')
title('K-means Clustering')
view(85,10)

for c=1:k  % PLOT CENTROID LOCATION
    %plot3(pc3(end-k+c,3),pc3(end-k+c,2),pc3(end-k+c,1),'+k','Markersize',50);
    plot3(pc3(end-k+c,3),pc3(end-k+c,2),pc3(end-k+c,1),'+','Markersize',30,'Color',cmap(c,:),'Linewidth',5);
end

% 
% % SAVE VIDEO Uncomment the lines below if you wish to save the video
% v = VideoWriter('Video_ClusterClouds.avi');
% v.FrameRate=5;
% open(v);
% EL=15;
% 
% for AZ=-180:10:180
%     % Set light and orientation angles
% %     AZ=Azymuth(t);
%     subplot(3,2,1)
%     view(AZ,EL)
%     
%     subplot(3,2,2)
%     view(AZ,EL)   
%     
%     pause(0.01)
%     
%         frame = getframe(gcf);
%         writeVideo(v,frame);
%         
% end
%     
% close(v)
%     
%     % VIDEO Here saves the current frame to the VIDEO
%     %     frame = getframe(gcf);
%     %     writeVideo(v,frame);
