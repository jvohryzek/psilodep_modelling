function [iFC, Var_Eig, Leading_Eig,iFCtril] = instFC( Phase_BOLD, numAreas, tMax,typeCalculation,typeDynamics )
        iFC = zeros(numAreas,numAreas,tMax);
        Isubdiag = find(tril(ones(numAreas),-1));
        
        %% for t=1 %%
    switch typeCalculation
        case 'default'   
            iFC(:,:,1) = cos(Phase_BOLD(:,1) - Phase_BOLD(:,1)');
            tmp = iFC(:,:,1);
            iFCtril(1,:) = tmp(Isubdiag);
            % Get the leading eigenvector
            [Vec1,Eig1] = eigs(iFC(:,:,1),1); % calculating all eigenvetors/values to obtain the variance explained (86 to make sure) BUT only 10 should do as most will be zeros
            % Calculate the variance explained by the leading eigenvector
            Var_Eig(1) = Eig1(1,1)/sum(abs(diag(Eig1)));
            % Make sure the leading eigenvector is the most negative WHY?
            if mean(Vec1>0)>.5 % if the number of regions more positive => make it negative
                Vec1 = -Vec1;
            elseif mean(Vec1>0) == .5 && sum(Vec1(Vec1>0))>-sum(Vec1(Vec1<0)) % if the number of regions between positive and negative is equal AND the sum of dominant weight is for the positive group => make it negative
                Vec1 = -Vec1;
            end
            Leading_Eig(:,1) = Vec1;
        
        case 'permute'
            % permuted case
            for n = 1:numAreas
                Phase_BOLD_permuted(n,:) = Phase_BOLD(n,randperm(231));
            end
                for n = 1:numAreas
                    for m = 1:numAreas
                        iFC(n,m,1) = cos(Phase_BOLD_permuted(n,1) - Phase_BOLD_permuted(m,1));
                    end
                end
            % Get the leading eigenvector
            [Vec1,Eig1] = eigs(iFC(:,:,1),2); % calculating all eigenvetors/values to obtain the variance explained (86 to make sure) BUT only 10 should do as most will be zeros
            Vec1 = Vec1(:,1);
            % Calculate the variance explained by the leading eigenvector
            Var_Eig(1) = Eig1(1,1)/sum(abs(diag(Eig1)));

            % Make sure the leading eigenvector is the most negative WHY?
                if mean(Vec1>0)>.5 % if the number of regions more positive => make it negative
                    Vec1 = -Vec1;
                elseif mean(Vec1>0) == .5 && sum(Vec1(Vec1>0))>-sum(Vec1(Vec1<0))% if the number of regions between positive and negative is equal AND the sum of dominant weight is for the positive group => make it negative
                    Vec1 = -Vec1;
                end
            Leading_Eig(:,1) = Vec1; 
    end
        %% for t=t+1 %%
        
        for t = 2:tMax

           switch typeCalculation
                case 'default'   
                    iFC(:,:,t) = cos(Phase_BOLD(:,t) - Phase_BOLD(:,t)');
                    tmp=iFC(:,:,t);
                    iFCtril(t,:) = tmp(Isubdiag);
                    % Get the leading eigenvector
                    [Vec1,Eig1] = eigs(iFC(:,:,t),1); % calculating all eigenvetors/values to obtain the variance explained (86 to make sure) BUT only 10 should do as most will be zeros
                    % calculating the previous timepoint
%                     [Vec1previous,~] = eigs(iFC(:,:,t-1),10); % previous timepoint
%                     Vec1previous = Vec1previous(:,1);
                    % Calculate the variance explained by the leading eigenvector
                    Var_Eig(t) = Eig1(1,1)/sum(abs(diag(Eig1)));
               switch typeDynamics
                   case 'halfSwitching'
                    % Make sure the leading eigenvector is the most negative WHY?
                        if mean(Vec1>0)>.5 % if the number of regions more positive => make it negative
                            Vec1 = -Vec1;
                        elseif mean(Vec1>0) == .5 && sum(Vec1(Vec1>0))>-sum(Vec1(Vec1<0))% if the number of regions between positive and negative is equal AND the sum of dominant weight is for the positive group => make it negative
                            Vec1 = -Vec1;
                        end
                    Leading_Eig(:,t) = Vec1;
                   case 'continuous'
                       % cosine distance => -1 opposite vectors, 1 exactly
                       % the same vectors, 0 orthogonal vectors
                       distV1vsV1previous = dot(Vec1,Leading_Eig(:,t-1))./(sqrt(sum(Vec1.^2)).*sqrt(sum(Leading_Eig(:,t-1).^2)));
                       %pdist([Vec1,Vec1previous]','cosine');
                       distmirrorV1vsV1previous = dot(-Vec1,Leading_Eig(:,t-1))./(sqrt(sum(Vec1.^2)).*sqrt(sum(Leading_Eig(:,t-1).^2)));
                       if distV1vsV1previous > distmirrorV1vsV1previous
                           Vec1 = Vec1;
                       elseif distmirrorV1vsV1previous > distV1vsV1previous
                           Vec1 = -Vec1;
                       end
                       Leading_Eig(:,t) = Vec1;
               end
                       
        
         case 'permute'
                    % permuted case
                    for n = 1:numAreas
                        Phase_BOLD_permuted(n,:) = Phase_BOLD(n,randperm(231));
                    end
                        for n = 1:numAreas
                            for m = 1:numAreas
                                iFC(n,m,t) = cos(Phase_BOLD_permuted(n,t) - Phase_BOLD_permuted(m,t));
                            end
                        end
                    % Get the leading eigenvector
                    [Vec1,Eig1] = eigs(iFC(:,:,t),10); % calculating all eigenvetors/values to obtain the variance explained (86 to make sure) BUT only 10 should do as most will be zeros
                    Vec1 = Vec1(:,1);
                    % Calculate the variance explained by the leading eigenvector
                    Var_Eig(t) = Eig1(1,1)/sum(abs(diag(Eig1)));
               switch typeDynamics
                   case 'halfSwitching'
                    % Make sure the leading eigenvector is the most negative WHY?
                        if mean(Vec1 > 0) > .5 % if the number of regions more positive => make it negative
                            Vec1 = -Vec1;
                        elseif mean(Vec1 > 0) == .5 && sum(Vec1(Vec1 > 0))>-sum(Vec1(Vec1 < 0))% if the number of regions between positive and negative is equal AND the sum of dominant weight is for the positive group => make it negative
                            Vec1 = -Vec1;
                        end
                    Leading_Eig(:,t) = Vec1;
                   case 'continuous'
                       distV1vsV1previous = pdist([Vec1,Vec1previous]','cosine');
                       distmirrorV1vsV1previous = pdist([Vec1,-Vec1previous]','cosine');
                       if distV1vsV1previous < distmirrorV1vsV1previous
                           Vec1 = Vec1;
                       elseif distV1vsV1previous > distmirrorV1vsV1previous
                           Vec1 = -Vec1;
                       end
                       Leading_Eig(:,t) = Vec1;
               end
   
                    Leading_Eig(:,t) = Vec1;               
            end 
        end

end

