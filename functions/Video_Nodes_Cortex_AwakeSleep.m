function Video_Nodes_Cortex_AwakeSleep %,IDX) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%    
%   Function to make a video of the Leading Eigenvector Dynamics 
%   in cortical space, where each sphere represents a brain area
%   colored according to the values of the eigenvector elements
%   at each frame. Links are between a proportion of the most connected
%   nodes.
%
%   Input Leading_Eig directly from Workspace or call inside a Function.
%   (or load Leading_Eig from a saved file)
%
%   Note: Leading_Eig had size (TxN) and the Leading Eigenvectors need to 
%   be corrected for symmetry. 
%
%   Writen by Joana Cabral
%   joana.cabral@psych.ox.ac.uk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SAVE VIDEO Uncomment the lines below if you wish to save the video
% v = VideoWriter('Video_LEiDA.avi');
% v.FrameRate=5;
% open(v);

load WakeSleep.mat Leading_Eig Time_all

Leading_Eig=Leading_Eig/max(abs(Leading_Eig(:)));
V1_awake=Leading_Eig((Time_all(2,:)==1),:);
V1_sleep=Leading_Eig((Time_all(2,:)==2),:);
clear Leading_Eig Time_all

% center origin
ori=[65 45.5 35];
Order=[1:2:90 90:-2:2];

load('/Users/joana/Documents/Work/LEiDA general/VideoLEiDA/aal_cog.txt','aal_cog')
MNI_coord=aal_cog(Order,:)/10;
clear aal_cog
scale=5.5;
[x,y,z] = sphere;
a=2;

figure

subplot(1,2,1)
hold on

% % PLOT CORTEX
cortex.pial=mapPial('MNI152_T1_2mm_brain_mask.nii');
cortex.val=0.2;
redux=1;
sregion=smooth3(cortex.pial);
psregion=patch(isosurface(sregion,cortex.val,'verbose'), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
reducepatch(psregion,redux,'verbose');
isonormals(sregion,psregion);
set(psregion,'FaceAlpha', 0.1); %transparency

% PLOT NODES
% Set here the first vector
Va=V1_awake(1,:);
for n=1:length(Va)
    if Va(n)>0
        ha(n)=surf(x*a+scale*MNI_coord(n,2)+ori(1), y*a+scale*MNI_coord(n,1)+ori(2),z*a+scale*MNI_coord(n,3)+ori(3),'FaceColor',[Va(n) 0 0],'EdgeColor','none','FaceAlpha',1);
  elseif Va(n)<0
        ha(n)=surf(x*a+scale*MNI_coord(n,2)+ori(1), y*a+scale*MNI_coord(n,1)+ori(2),z*a+scale*MNI_coord(n,3)+ori(3),'FaceColor',[0 0 abs(Va(n))],'EdgeColor','none','FaceAlpha',.5);
    end
end

axis off;
axis equal
material dull; lighting phong;
view([-90 90]) % top
l=camlight;

subplot(1,2,2)
hold on

% % PLOT CORTEX
cortex.pial=mapPial('MNI152_T1_2mm_brain_mask.nii');
cortex.val=0.2;
redux=1;
sregion=smooth3(cortex.pial);
psregion=patch(isosurface(sregion,cortex.val,'verbose'), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
reducepatch(psregion,redux,'verbose');
isonormals(sregion,psregion);
set(psregion,'FaceAlpha', 0.1); %transparency

% PLOT NODES
% Set here the first vector
Vs=V1_sleep(1,:);

for n=1:length(Vs)
    if Vs(n)>0
        hs(n)=surf(x*a+scale*MNI_coord(n,2)+ori(1), y*a+scale*MNI_coord(n,1)+ori(2),z*a+scale*MNI_coord(n,3)+ori(3),'FaceColor',[Vs(n) 0 0],'EdgeColor','none','FaceAlpha',1);
  elseif Vs(n)<0
        hs(n)=surf(x*a+scale*MNI_coord(n,2)+ori(1), y*a+scale*MNI_coord(n,1)+ori(2),z*a+scale*MNI_coord(n,3)+ori(3),'FaceColor',[0 0 abs(Vs(n))],'EdgeColor','none','FaceAlpha',.5);
    end
end

axis off;
axis equal
material dull; lighting phong;
view([-90 90]) % top
l=camlight;

%   EL=-10; % the desired elevation of the camera
%   Azymuth=-180:+1:+180; % To rotate the azymuth around
  
  %cmap=[0 0 1; 0 1 1 ; 0 1 0 ; 1 1 0;  1 0 0; 1 0 1; .7 .7 .7 ; 1 0.5 0 ];

for t=1:length(V1_awake)

    subplot(1,2,1)
    hold on
    Va=V1_awake(t,:);
    
    %Update here the colors of the spheres
    for n=1:length(Va)
        if Va(n)>0
            set(ha(n),'FaceColor',[Va(n) 0 0]);
        elseif Va(n)<0
            set(ha(n),'FaceColor',[0 0 abs(Va(n))]);
        end
    end
    
    cha=[];    
    n_strong=find(Va>0.5);
    if numel(n_strong)>1
        u=1;

        for a=1:numel(n_strong)
            n=n_strong(a);
            for b=1:a
                p=n_strong(b);
                c1=[scale*MNI_coord(n,2)+ori(1) scale*MNI_coord(n,1)+ori(2) scale*MNI_coord(n,3)+ori(3)];
                c2=[scale*MNI_coord(p,2)+ori(1) scale*MNI_coord(p,1)+ori(2) scale*MNI_coord(p,3)+ori(3)];
                cha(u)=plot3([c1(1) c2(1)],[c1(2) c2(2)],[c1(3) c2(3)],'Color','r'); %cmap(IDX(t),:));
                u=u+1;
            end
        end
    end
    
    chs=[];    

    subplot(1,2,2)
    hold on
    Vs=V1_sleep(t,:);
    
    %Update here the colors of the spheres
    for n=1:length(Vs)
        if Vs(n)>0
            set(hs(n),'FaceColor',[Vs(n) 0 0]);
        elseif Vs(n)<0
            set(hs(n),'FaceColor',[0 0 abs(Vs(n))]);
        end
    end
    
  
    n_strong=find(Vs>0.65);
    if numel(n_strong)>1
        u=1;

        for a=1:numel(n_strong)
            n=n_strong(a);
            for b=1:a
                p=n_strong(b);
                c1=[scale*MNI_coord(n,2)+ori(1) scale*MNI_coord(n,1)+ori(2) scale*MNI_coord(n,3)+ori(3)];
                c2=[scale*MNI_coord(p,2)+ori(1) scale*MNI_coord(p,1)+ori(2) scale*MNI_coord(p,3)+ori(3)];
                chs(u)=plot3([c1(1) c2(1)],[c1(2) c2(2)],[c1(3) c2(3)],'Color','r'); %cmap(IDX(t),:));
                u=u+1;
            end
        end
    end
   
    % Set light and orientation angles
%     AZ=Azymuth(t);
%     view(AZ,EL)
%     camorbit(10,0)
%     camlight(l)
    
    % VIDEO Here saves the current frame to the VIDEO
    %     frame = getframe(gcf);
    %     writeVideo(v,frame);
    
    % Add a pause to see the video on screen
    pause(0.1)
    
    % Delete all links in each frame
    delete(cha)
    delete(chs)
end

% SAVE VIDEO uncomment to close and finish the video
% This line finishes the video
% close(v)

end

% Function needed to plot the cortex

function pial=mapPial(region)

VG=spm_vol(region(1,:));
pial=zeros(VG.dim(1:3));
for i=1:VG.dim(3)
    pial(:,:,i) = spm_slice_vol(VG,spm_matrix([0 0 i]),VG.dim(1:2),1);
end

end