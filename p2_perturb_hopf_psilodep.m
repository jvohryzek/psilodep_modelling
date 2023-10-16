clear all;

load  empiricalLEiDA.mat;
load  optimizedhopfawake.mat

P1emp=mean(P1emp);
P2emp=mean(P2emp);

PERTURB=0:0.01:0.1;  %% This is synchronisation protocol
%PERTURB=0:-0.025:-0.4; %% This is for the noise Protocol
we=0.09;
ITER=30;
N=90;
a=zeros(N,2);
C=squeeze(Coptim(find(abs(WE-we)<0.0001),:,:));

TSmax=1000;
NSUB=n_Subjects;
TR=2.08;  % Repetition Time (seconds)
NumClusters=Number_Clusters;

%%%%%%%%%%%%%%%%%%

omega = repmat(2*pi*f_diff',1,2); omega(:,1) = -omega(:,1);

dt=0.1*TR/2;
Tmax=NSUB*TSmax;
sig=0.02;
dsig = sqrt(dt)*sig; % to avoid sqrt(dt) at each time step

%%%%%%%%%%%%
%% Optimize
%%
wC = we*C;
sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj

for node=1:N/2
    node
    iwe=1;
    for perturb=PERTURB
        a=zeros(N,2);
        a(node,:)=a(node,:)+perturb;
        a(N+1-node,:)=a(N+1-node,:)+perturb;
        for iter=1:ITER
            xs=zeros(Tmax,N);
            z = 0.1*ones(N,2); % --> x = z(:,1), y = z(:,2)
            nn=0;
            % discard first 3000 time steps
            for t=0:dt:3000
                suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
                zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
                z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
            end
            % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
            for t=0:dt:((Tmax-1)*TR)
                suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
                zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
                z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
                if abs(mod(t,TR))<0.01
                    nn=nn+1;
                    xs(nn,:)=z(:,1)';
                end
            end
            
            %%%% KL dist between PTR2emp
            
            [PTRsim,Pstates]=LEiDA_fix_cluster(xs',NumClusters,Vemp,TR);
            KLps_p(iter)=0.5*(sum(Pstates.*log(Pstates./P2emp))+sum(P2emp.*log(P2emp./Pstates)));
        end
        KLpstatessleep_perturbed(node,iwe)=mean(KLps_p)
        iwe=iwe+1;
    end
end

imagesc(KLpstatessleep_perturbed);

