function [OP Synchro Metasta] = KuramotoOP( Phase_BOLD, N_areas)
        
% Calculate the Kuramoto Order Parameter

OP = abs(sum(exp(1i * Phase_BOLD))/N_areas);
Metasta = std(OP);
Synchro = mean(OP);

end

