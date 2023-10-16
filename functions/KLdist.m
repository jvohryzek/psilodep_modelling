function dist=KLdist(TRANS_1,TRANS_2)
ns=size(TRANS_1,1);
EMIS_1  = eye(ns);
EMIS_2  = eye(ns);

% check row full of zeros
%--------------------------------------
flag_zeros = 0;
S1 = sum(TRANS_1,2);
S2 = sum(TRANS_2,2);
if sum(S1)<ns
    flag_zeros = 1; 
elseif sum(S2)<ns
    flag_zeros = 2; 
end

if flag_zeros == 1
        ind=find(S1==0);
        j = setxor(1:ns,ind);
        j = j(1);
        for i=ind
            TRANS_1(i,i)=0;
            TRANS_1(i,j)=1;
        end
elseif flag_zeros == 2
        ind=find(S2==0);
        j = setxor(1:ns,ind);
        j = j(1);
        for i=ind
            TRANS_2(i,i)=0;
            TRANS_2(i,j)=1;
        end        
end
%--------------------------------------

% Check ergodicity:
if any(TRANS_1(:)==0)
    [i,j] = find(TRANS_1==0);
    rows  = unique(i);
    for k=1:length(rows)
        cols = j(i==rows(k));
        q = length(cols);
        TR = TRANS_1(k,:);
        [~,jmax] = max(TR);
        TR(cols) = eps;
        TR(jmax) = TR(jmax) - q*eps;
        TRANS_1(rows(k),:) = TR;
    end
end

if any(TRANS_2(:)==0)
    [i,j] = find(TRANS_2==0);
    rows  = unique(i);
    for k=1:length(rows)
        cols = j(i==rows(k));
        q = length(cols);
        TR = TRANS_2(k,:);
        [~,jmax] = max(TR);
        TR(cols) = eps;
        TR(jmax) = TR(jmax) - q*eps;
        TRANS_2(rows(k),:) = TR;
    end
end
    

% Here starts the comparison:

% A long observation generated by model 1
T = 1000000; 
O1 = hmmgenerate(T,TRANS_1,EMIS_1);

% A long observation generated by model 2
T = 1000000; 
O2 = hmmgenerate(T,TRANS_2,EMIS_2);

% get log-likelihoods ( logP ):

[~,logP11] = hmmdecode(O1,TRANS_1,EMIS_1);
[~,logP12] = hmmdecode(O1,TRANS_2,EMIS_2);
[~,logP21] = hmmdecode(O2,TRANS_1,EMIS_1);
[~,logP22] = hmmdecode(O2,TRANS_2,EMIS_2);

% Divergences:

D12 = ( logP11 - logP12 );

D21 = ( logP22 - logP21 );

% symmetrize:

dist = 1/2*( D12 + D21 )/T;
