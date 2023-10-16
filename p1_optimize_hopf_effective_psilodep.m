clear all;

load empiricalLEiDA.mat;

P1emp=mean(P1emp);
P2emp=mean(P2emp);

load sc90.mat;
C=sc90;
C=C/max(max(C))*0.2;

load data_Awake.mat;

TSmax=240;
NSUB=18;
TR=2.08;  % Repetition Time (seconds)
NumClusters=Number_Clusters;

delt = TR;            % sampling interval
k=2;                  % 2nd order butterworth filter
fnq=1/(2*delt);
flp = .04;           % lowpass frequency of filter
fhi = fnq-0.001;           % highpass
Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt,afilt]=butter(k,Wn);   % construct the filter

flp = .04;           % lowpass frequency of filter
fhi = .07;           % highpass
Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt2,afilt2]=butter(k,Wn);   % construct the filter
clear fnq flp fhi Wn k

for nsub=1:NSUB
    [N, Tmax0]=size(X{1,nsub});
    if Tmax0<TSmax
        X{1,nsub}=[];
    end
end
X=X(~cellfun('isempty',X));
n_Subjects=size(X,2);

%%%%%%%%%%%%%%
kk=1;
insub=1;
TSmax=1000;

N=90;
Isubdiag = find(tril(ones(N),-1));

Tmaxtotal=0;
for nsub=1:n_Subjects
    [N, Tmax0]=size(X{1,nsub});
    Tmax=min(TSmax,Tmax0);
    Tmaxtotal=Tmaxtotal+Tmax;
    signaldata = X{1,nsub};
    signaldata=signaldata(:,1:Tmax);
    Phase_BOLD_data=zeros(N,Tmax);
    timeseriedata=zeros(N,Tmax);
    for seed=1:N
        x=demean(detrend(signaldata(seed,:)));
        x(find(x>3*std(x)))=3*std(x);
        x(find(x<-3*std(x)))=-3*std(x);
        timeseriedata(seed,:) = filtfilt(bfilt2,afilt2,x);    % zero phase filter the data
        Phase_BOLD_data(seed,:) = angle(hilbert(timeseriedata(seed,:)));
    end
    T=10:Tmax-10;
    for t=T
        kudata=sum(complex(cos(Phase_BOLD_data(:,t)),sin(Phase_BOLD_data(:,t))))/N;
        syncdata(t-9)=abs(kudata);
        for i=1:N
            for j=1:i-1
                patt(i,j)=cos(adif(Phase_BOLD_data(i,t),Phase_BOLD_data(j,t)));
            end
        end
        pattern(t-9,:)=patt(Isubdiag);
    end
    metastabilitydata2(nsub)=std(syncdata);
    
    npattmax=size(pattern,1);
    for t=1:npattmax-2
        p1=mean(pattern(t:t+2,:));
        for t2=t+1:npattmax-2
            p2=mean(pattern(t2:t2+2,:));
            phfcddata(kk)=dot(p1,p2)/norm(p1)/norm(p2);
            kk=kk+1;
        end
    end
    
    for t=1:Tmax
        for n=1:N
            for p=1:N
                iFC(t,n,p)=cos(Phase_BOLD_data(n,t)-Phase_BOLD_data(p,t));
            end
        end
    end
    FCphasesemp2(nsub,:,:)=squeeze(mean(iFC));
end
FCphasesemp=squeeze(mean(FCphasesemp2));
metastabilitydata=mean(metastabilitydata2);


for nsub=1:n_Subjects
    clear PowSpect PowSpect2;
    [N, Tmax0]=size(X{1,nsub});
    Isubdiag = find(tril(ones(N),-1));
    Tmax=min(TSmax,Tmax0);
    TT=Tmax;
    Ts = TT*TR;
    freq = (0:TT/2-1)/Ts;
    signaldata = X{1,nsub};
    signaldata=signaldata(:,1:Tmax);
    FCemp2(nsub,:,:)=corrcoef(signaldata');
    
    %%%%
    
    [aux minfreq]=min(abs(freq-0.04));
    [aux maxfreq]=min(abs(freq-0.07));
    nfreqs=length(freq);
    
    
    for seed=1:N
        x=detrend(demean(signaldata(seed,:)));
        ts =zscore(filtfilt(bfilt2,afilt2,x));
        pw = abs(fft(ts));
        PowSpect(:,seed,insub) = pw(1:floor(TT/2)).^2/(TT/TR);
        ts2 =zscore(filtfilt(bfilt,afilt,x));
        pw2 = abs(fft(ts2));
        PowSpect2(:,seed,insub) = pw2(1:floor(TT/2)).^2/(TT/TR);
    end
    insub=insub+1;
end

Power_Areas=mean(PowSpect,3);
Power_Areas2=mean(PowSpect2,3);
for seed=1:N
    Power_Areas(:,seed)=gaussfilt(freq,Power_Areas(:,seed)',0.01);
    Power_Areas2(:,seed)=gaussfilt(freq,Power_Areas2(:,seed)',0.01);
    vsig(seed)=sum(Power_Areas2(minfreq:maxfreq,seed))/sum(Power_Areas2(:,seed));
end

vmax=max(vsig);
vmin=min(vsig);

[maxpowdata,index]=max(Power_Areas);
f_diff = freq(index);
FCemp=squeeze(mean(FCemp2));

clear PowSpect PowSpect2 Power_Areas Power_Areas2;

%%%%%%%%%%%%%%%%%%

omega = repmat(2*pi*f_diff',1,2); omega(:,1) = -omega(:,1);

dt=0.1*TR/2;
Tmax=TSmax*n_Subjects;
sig=0.02;
dsig = sqrt(dt)*sig; % to avoid sqrt(dt) at each time step

%%%%%%%%%%%%
%% Optimize
%%
iwe=1;
WE=0:0.01:1;%0:0.005:0.15;
a=zeros(N,2);

NWE=length(WE);
PTRsimul=zeros(NWE,NumClusters,NumClusters);
Pstatessimul=zeros(NWE,NumClusters);

for we=WE
    minm=100;
    Cnew=C;
    for iter=1:150
        wC = we*Cnew;
        sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj
        xs=zeros(Tmax,N);
        %number of iterations, 100 willk�hrlich, weil reicht in diesem Fall
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
        
        %%%%
        BOLD=xs';
        signal_filt=zeros(N,nn);
        Phase_BOLD=zeros(N,nn);
        for seed=1:N
            BOLD(seed,:)=demean(detrend(BOLD(seed,:)));
            signal_filt(seed,:) =filtfilt(bfilt2,afilt2,BOLD(seed,:));
            Phase_BOLD(seed,:) = angle(hilbert(signal_filt(seed,:)));
        end
        
        for t=1:nn
            for n=1:N
                for p=1:N
                    iFC(t,n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
                end
            end
        end
        FCphases=squeeze(mean(iFC));
        
        for i=1:N
            for j=i+1:N
               % if (C(i,j)>0 || j==N-i+1)
                    Cnew(i,j)=Cnew(i,j)+0.01*(FCphasesemp(i,j)-FCphases(i,j));
                    if Cnew(i,j)<0
                        Cnew(i,j)=0;
                    end
                    Cnew(j,i)=Cnew(i,j);
               % end
            end
        end
        
        Cnew=Cnew/max(max(Cnew))*0.2;
        
        D = abs(FCphasesemp-FCphases).^2;
        MSE = sum(D(:))/numel(FCphases);
        if MSE<0.01
            break;
        end
        
        %%%%

    end
        
    Coptim(iwe,:,:)=Cnew;

    %%%%%%%%%%%%%%
    %%% Final simul
    
    xs=zeros(Tmax,N);
    %number of iterations, 100 willk�hrlich, weil reicht in diesem Fall
    z = 0.1*ones(N,2); % --> x = z(:,1), y = z(:,2)
    nn=0;
    % discard first 3000 time steps
    for t=0:dt:3000
        suma = wC*z - sumC.*z; % where SC enters: sum(Cij*xi) - sum(Cij)*xj
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
    
    FC_simul=corrcoef(xs(1:nn,:));
    cc=corrcoef(atanh(FCemp(Isubdiag)),atanh(FC_simul(Isubdiag)));
    fitt(iwe)=cc(2);
    
    %%%%% Meta & FCD
    BOLD=xs';
    Phase_BOLD=zeros(N,nn);
    signal_filt=zeros(N,nn);
    for seed=1:N
        BOLD(seed,:)=demean(detrend(BOLD(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt2,afilt2,BOLD(seed,:));
        Phase_BOLD(seed,:) = angle(hilbert(signal_filt(seed,:)));
    end
    T=10:Tmax-10;
    for t=T
        ku=sum(complex(cos(Phase_BOLD(:,t)),sin(Phase_BOLD(:,t))))/N;
        sync(t-9)=abs(ku);
        for i=1:N
            for j=1:i-1
                patt(i,j)=cos(adif(Phase_BOLD(i,t),Phase_BOLD(j,t)));
            end
        end
        pattern(t-9,:)=patt(Isubdiag);
    end
    metastability(iwe)=abs(metastabilitydata-std(sync));
    
    kk=1;
    npattmax=size(pattern,1);
    for t=1:npattmax-2
        p1=mean(pattern(t:t+2,:));
        for t2=t+1:npattmax-2
            p2=mean(pattern(t2:t2+2,:));
            phfcd(kk)=dot(p1,p2)/norm(p1)/norm(p2);
            kk=kk+1;
        end
    end
    
    [H,P,ksdist(iwe)]=kstest2(phfcd,phfcddata);
    
    %%%% KL dist between PTR2emp
    
    [PTRsim,Pstates]=LEiDA_fix_cluster(xs',NumClusters,Vemp,TR);
    
    klpstatessleep(iwe)=0.5*(sum(Pstates.*log(Pstates./P2emp))+sum(P2emp.*log(P2emp./Pstates)));
    klpstatesawake(iwe)=0.5*(sum(Pstates.*log(Pstates./P1emp))+sum(P1emp.*log(P1emp./Pstates)));

    kldistsleep(iwe)=KLdist(PTR2emp,PTRsim);
    kldistawake(iwe)=KLdist(PTR1emp,PTRsim);
    entropydistsleep(iwe)=EntropyMarkov(PTR2emp,PTRsim);
    entropydistawake(iwe)=EntropyMarkov(PTR1emp,PTRsim);
    
    PTRsimul(iwe,:,:)=PTRsim;
    Pstatessimul(iwe,:)=Pstates;
    
    iwe=iwe+1;
    
    ksdist
    klpstatesawake

end

save optimizedhopfawake.mat WE PTRsimul Pstatessimul metastability ksdist klpstatessleep klpstatesawake kldistsleep kldistawake entropydistawake entropydistsleep fitt Coptim n_Subjects f_diff;

figure
plot(WE,fitt,'b');
hold on;
plot(WE,kldistawake,'k');
plot(WE,entropydistawake,'k');    
figure
plot(WE,metastability,'r');
hold on;
plot(WE,ksdist,'c');
plot(WE,klpstatesawake,'k');
plot(WE,klpstatessleep,'b');

% 
% figure
% plot(WE,fitt,'b');
% figure
% plot(WE,kldistsleep,'r');
% hold on;
% plot(WE,kldistawake,'k');
% figure
% plot(WE,entropydistsleep,'r');
% hold on;
% plot(WE,entropydistawake,'k');    
% figure
% plot(WE,klpstatessleep,'r');
% hold on;
% plot(WE,klpstatesawake,'k');   
        