function [nmiOpt, purityOpt, OAOpt, AAOpt, kappaOpt, timeIdx]= measure_performance(Clusterings, Y)
 
if isstruct(Clusterings)
    % Multiple clusterings. We store the performance of the best one.

    % Calculate the NMI and purity between best clustering and GT  
    nmiTemp = zeros(length(Clusterings.K),1);
    purityTemp = zeros(length(Clusterings.K),1);
    OATemp = zeros(length(Clusterings.K),1);
    kappaTemp = zeros(length(Clusterings.K),1);
    AATemp = zeros(length(Clusterings.K),1);
    for t = 1:length(purityTemp)
        % Restrict analysis to points with labels>1. This corresponds to
        % labeled ground truth points.
        nmiTemp(t) = nmi(double(Clusterings.Labels(Y>1,t)), double(Y(Y>1)));
        purityTemp(t) = purity(double(Clusterings.Labels(Y>1,t)), double(Y(Y>1)));

        % Unlabeled points are excluded in this analysis already.
        [ OATemp(t), kappaTemp(t), AATemp(t)] = calcAccuracy(double(Y), double(Clusterings.Labels(:,t)));

    end

    [OAOpt, timeIdx] = max(OATemp);
    AAOpt = AATemp(timeIdx);
    kappaOpt = kappaTemp(timeIdx);
    nmiOpt = nmiTemp(timeIdx);
    purityOpt = purityTemp(timeIdx);
 

elseif sum(size(Clusterings) == size(Y)) == 2
    
    timeIdx = 1;
    [ OAOpt, kappaOpt, AAOpt] = calcAccuracy(Y, Clusterings);
    nmiOpt = nmi(Clusterings(Y>1), Y(Y>1));
    purityOpt = purity(Clusterings(Y>1), Y(Y>1));
     
end