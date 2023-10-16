function [hAdd2] = plotKSdistStats(statsKS)
% Calculates the KS statistics plot
% hAdd2 => KSttest pvalues as a function of bins
% hAdd3 => KSttest values as a function of bins

%%%%%%%%%%%%% plotting p and ktest for the hist. difference %%%%%%%%%%%


% pvalue
hAdd2 = figure; opt.FontSize = 18;
subplot(2,1,1)
plot(10:100:10000,statsKS.preal,'LineWidth',2)
hold on
plot(10:100:10000,statsKS.prealrem,'LineWidth',2)
hold on
plot(10:100:10000,statsKS.prealrembfvsafr,'.-','LineWidth',2)
hold on
plot(10:100:10000,statsKS.prealrembfvsafnr,'.-','LineWidth',2)
hold on
plot(10:100:10000,statsKS.prand,'--','LineWidth',2)
hold on
plot(10:100:10000,statsKS.prandshifted,'--','LineWidth',2)
hold on
plot(10:100:10000,statsKS.prealrembfvsbfnr,'-','LineWidth',2)
title('two-sample KS test with varying bin size'); xlabel('Bin Size'); ylabel('p-value');
set(gca,'FontSize',24); legend({'before/after','after responders/non-responders',...
    'before all vs. after responders','before all vs. after non-responders',...
    'normal','a*normal+b','before responders/non-responders'});

% kS-distance
subplot(2,1,2)
plot(10:100:10000,statsKS.ktestreal,'LineWidth',2)
hold on
plot(10:100:10000,statsKS.ktestrealrem,'LineWidth',2)
hold on
plot(10:100:10000,statsKS.prealrembfvsafr,'.-','LineWidth',2)
hold on
plot(10:100:10000,statsKS.prealrembfvsafnr,'.-','LineWidth',2)
hold on
plot(10:100:10000,statsKS.ktestrand,'--','LineWidth',2)
hold on
plot(10:100:10000,statsKS.ktestrandshifted,'--','LineWidth',2)
hold on
plot(10:100:10000,statsKS.prealrembfvsbfnr,'-','LineWidth',2)
title('two-sample KS test with varying bin size'); xlabel('Bin Size'); ylabel('KS-distance')
set(gca,'FontSize',24); legend({'before/after','after responders/non-responders',...
    'before all vs. after responders','before all vs. after non-responders',...
    'normal','a*normal+b','before responders/non-responders'});
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
hAdd2 = fancy_figure(hAdd2, opt);

end

