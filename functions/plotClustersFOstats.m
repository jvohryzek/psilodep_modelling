function [g] = plotClustersFOstats(Cluster,typeTest,typeMatrix)
%This function plots all the brain clusters fro the varying cluster value
%sorted based on the cluster's largest fractional occupancy
%
switch typeTest
    case 'ttest'
    g = figure;opt.FontSize = 18;
    subplot(2,6,1);imagesc(tril(Cluster.hRest,0));title('After/Before'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24)
    subplot(2,6,2);imagesc(tril(Cluster.hRestNoResp,0));title('After/Before NoResp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24)
    subplot(2,6,3);imagesc(tril(Cluster.hRestResp,0));title('After/Before Resp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24)
    subplot(2,6,4);imagesc(tril(Cluster.hRestdiff,0));title(' (After - Before)Resp/NoResp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24)
    subplot(2,6,5);imagesc(tril(Cluster.hRestAfRestNoRest,0));title(' After Resp/NoResp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24) 
    subplot(2,6,6);imagesc(tril(Cluster.hRestBfAllAfRest,0));title(' Before All/ After Resp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24) 

    subplot(2,6,7);imagesc(tril(Cluster.pRest,0));title('After - Before'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24)
    subplot(2,6,8);imagesc(tril(Cluster.pRestNoResp,0));title('After - Before NoResp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24)
    subplot(2,6,9);imagesc(tril(Cluster.pRestResp,0));title('After - Before Resp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24)
    subplot(2,6,10);imagesc(tril(Cluster.pRestdiff,0));title(' (After - Before)Resp/NoResp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24) 
    subplot(2,6,11);imagesc(tril(Cluster.pRestAfRestNoRest,0));title(' After Resp/NoResp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24) 
    subplot(2,6,12);imagesc(tril(Cluster.pRestBfAllAfRest,0));title(' Before All/ After Resp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24) 
    suptitle('Ttest')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    g = fancy_figure(g, opt);
    case 'ttestperm'
        switch typeMatrix
            case 'upperTri'
                g = figure;opt.FontSize = 18;
                subplot(2,6,1);imagesc(tril(Cluster.hRest,0));title('After/Before'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,2);imagesc(tril(Cluster.hRestNoResp,0));title('After/Before NoResp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,3);imagesc(tril(Cluster.hRestResp,0));title('After/Before Resp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,4);imagesc(tril(Cluster.hRestdiff,0));title(' (After - Before)Resp/NoResp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,5);imagesc(tril(Cluster.hRestAfRestNoRest,0));title(' After Resp/NoResp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,6);imagesc(tril(Cluster.hRestBfAllAfRest,0));title(' Before All/ After Resp'); xlabel('cluster');ylabel('cluster');

                subplot(2,6,7);imagesc(tril(Cluster.pRest,0));title('After - Before'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,8);imagesc(tril(Cluster.pRestNoResp,0));title('After - Before NoResp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,9);imagesc(tril(Cluster.pRestResp,0));title('After - Before Resp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,10);imagesc(tril(Cluster.pRestdiff,0));title(' (After - Before)Resp/NoResp'); xlabel('cluster');ylabel('cluster'); 
                subplot(2,6,11);imagesc(tril(Cluster.pRestAfRestNoRest,0));title(' After Resp/NoResp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,12);imagesc(tril(Cluster.pRestBfAllAfRest,0));title(' Before All/ After Resp'); xlabel('cluster');ylabel('cluster');

            case 'wholeMatrix'
                g = figure;opt.FontSize = 18;
                subplot(2,6,1);imagesc(Cluster.hRest);title('After/Before'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,2);imagesc(Cluster.hRestNoResp);title('After/Before NoResp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,3);imagesc(Cluster.hRestResp);title('After/Before Resp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,4);imagesc(Cluster.hRestdiff);title(' (After - Before)Resp/NoResp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,5);imagesc(Cluster.hRestAfRestNoRest);title(' After Resp/NoResp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,6);imagesc(Cluster.hRestBfAllAfRest);title(' Before All/ After Resp'); xlabel('cluster');ylabel('cluster');

                subplot(2,6,7);imagesc(Cluster.pRest);title('After - Before'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,8);imagesc(Cluster.pRestNoResp);title('After - Before NoResp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,9);imagesc(Cluster.pRestResp);title('After - Before Resp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,10);imagesc(Cluster.pRestdiff);title(' (After - Before)Resp/NoResp'); xlabel('cluster');ylabel('cluster'); 
                subplot(2,6,11);imagesc(Cluster.pRestAfRestNoRest);title(' After Resp/NoResp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,12);imagesc(Cluster.pRestBfAllAfRest);title(' Before All/ After Resp'); xlabel('cluster');ylabel('cluster');
             end
    suptitle('Permuted ttest')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    g = fancy_figure(g, opt);
    case 'ranksumperm'
        switch typeMatrix
            case 'upperTri'
                g = figure;opt.FontSize = 18;
                subplot(2,6,1);imagesc(tril(Cluster.hNPRest,0));title('After/Before'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,2);imagesc(tril(Cluster.hNPRestNoResp,0));title('After/Before NoResp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,3);imagesc(tril(Cluster.hNPRestResp,0));title('After/Before Resp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,4);imagesc(tril(Cluster.hNPRestdiff,0));title(' (After - Before)Resp/NoResp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,5);imagesc(tril(Cluster.hNPRestAfRestNoRest,0));title(' After Resp/NoResp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,6);imagesc(tril(Cluster.hNPRestBfAllAfRest,0));title(' Before All/ After Resp'); xlabel('cluster');ylabel('cluster');

                subplot(2,6,7);imagesc(tril(Cluster.pNPRest,0));title('After - Before'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,8);imagesc(tril(Cluster.pNPRestNoResp,0));title('After - Before NoResp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,9);imagesc(tril(Cluster.pNPRestResp,0));title('After - Before Resp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,10);imagesc(tril(Cluster.pNPRestdiff,0));title(' (After - Before)Resp/NoResp'); xlabel('cluster');ylabel('cluster'); 
                subplot(2,6,11);imagesc(tril(Cluster.pNPRestAfRestNoRest,0));title(' After Resp/NoResp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,12);imagesc(tril(Cluster.pNPRestBfAllAfRest,0));title(' Before All/ After Resp'); xlabel('cluster');ylabel('cluster');

                case 'wholeMatrix'
                g = figure;opt.FontSize = 18;
                subplot(2,6,1);imagesc(Cluster.hNPRest);title('After/Before'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,2);imagesc(Cluster.hNPRestNoResp);title('After/Before NoResp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,3);imagesc(Cluster.hNPRestResp);title('After/Before Resp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,4);imagesc(Cluster.hNPRestdiff);title(' (After - Before)Resp/NoResp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,5);imagesc(Cluster.hNPRestAfRestNoRest);title(' After Resp/NoResp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,6);imagesc(Cluster.hNPRestBfAllAfRest);title(' Before All/ After Resp'); xlabel('cluster');ylabel('cluster');

                subplot(2,6,7);imagesc(Cluster.pNPRest);title('After - Before'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,8);imagesc(Cluster.pNPRestNoResp);title('After - Before NoResp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,9);imagesc(Cluster.pNPRestResp);title('After - Before Resp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,10);imagesc(Cluster.pNPRestdiff);title(' (After - Before)Resp/NoResp'); xlabel('cluster');ylabel('cluster'); 
                subplot(2,6,11);imagesc(Cluster.pNPRestAfRestNoRest);title(' After Resp/NoResp'); xlabel('cluster');ylabel('cluster');
                subplot(2,6,12);imagesc(Cluster.pNPRestBfAllAfRest);title(' Before All/ After Resp'); xlabel('cluster');ylabel('cluster');
        end
        suptitle('Permuted ranksum')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        g = fancy_figure(g, opt);



    case 'ranksum'
    g = figure;opt.FontSize = 18;
    subplot(2,6,1);imagesc(tril(Cluster.hRestWsignrank,0));title('After/Before'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24)
    subplot(2,6,2);imagesc(tril(Cluster.hRestNoRespWsignrank,0));title('After/Before NoResp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24)
    subplot(2,6,3);imagesc(tril(Cluster.hRestRespWsignrank,0));title('After/Before Resp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24)
    subplot(2,6,4);imagesc(tril(Cluster.hRestdiffWranksum,0));title(' (After - Before)Resp/NoResp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24)
    subplot(2,6,5);imagesc(tril(Cluster.hRestAfRestNoRespWranksum,0));title(' After Resp/NoResp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24) 
    subplot(2,6,6);imagesc(tril(Cluster.hRestBfAllAfRestWranksum,0));title(' Before All/ After Resp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24) 

    subplot(2,6,7);imagesc(tril(Cluster.pRestWsignrank,0));title('After - Before'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24)
    subplot(2,6,8);imagesc(tril(Cluster.pRestNoRespWsignrank,0));title('After - Before NoResp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24)
    subplot(2,6,9);imagesc(tril(Cluster.pRestRespWsignrank,0));title('After - Before Resp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24)
    subplot(2,6,10);imagesc(tril(Cluster.pRestdiffWranksum,0));title(' (After - Before)Resp/NoResp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24) 
    subplot(2,6,11);imagesc(tril(Cluster.pRestAfRestNoRespWranksum,0));title(' After Resp/NoResp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24) 
    subplot(2,6,12);imagesc(tril(Cluster.pRestBfAllAfRestWranksum,0));title(' Before All/ After Resp'); xlabel('cluster');ylabel('cluster');set(gca,'FontSize',24) 

    suptitle('ranksum')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    g = fancy_figure(g, opt);
end
end

