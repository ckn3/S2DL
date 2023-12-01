function [ OA, kappa, tIdx, OATemp, kappaTemp, AA] = calcPerformance(Y, C, ignore1flag)
%{ 
Calculates statistics on clusterings
Inputs: Y (GT Labels), C (Estimated Clustering), and ignore1flag (1 if Y=1
        class is ignored in performance calculations).
Outputs: OA (Overall Accuracy), kappa (Cohen's kappa), tIdx (optimal
         clustering index), OATemp (OA values for each tIdx) and kappaTemp 
         (kappa values for each tIdx)
C may be one of the following formats: 
    - Structure with "Labels" field that is an nxM array with a clustering
    of X in each column. 
    - an nx1 clustering of X.
Copyright: Sam L. Polk (2022).
%}
if isstruct(C)

    numC = size(C.Labels,2);
    OATemp = zeros(numC,1);
    kappaTemp = zeros(numC,1);
    AATemp = zeros(numC,1);
    for i = 1:numC
        [ OATemp(i), kappaTemp(i), AATemp(i)] = evaluatePerformances(Y, C.Labels(:,i), ignore1flag);
    end
    [~, tIdx] = max(OATemp+kappaTemp+AATemp);
    OA = OATemp(tIdx);
    kappa = kappaTemp(tIdx);
    AA = AATemp(tIdx);

else
    % Single clustering (no need to run over multiple clusterings
    [ OA, kappa, AA] = evaluatePerformances(Y, C, ignore1flag);
    tIdx =1;
    OATemp = 0;
    kappaTemp = 0;
end

end

function [OA, kappa, AA] = evaluatePerformances(Y,C, ignore1flag)
    if ignore1flag
        % If true, we restric performance evaluation to unlabeled points (those
        % marked as index 1).
    
        % Perform hungarian algorithm to align clustering labels
        CNew = C(Y>1);
    
        missingk = setdiff(1:max(CNew), unique(CNew)');
        if length(missingk) == 1
            CNew(CNew>=missingk) = CNew(CNew>=missingk) - 1;
        else
            Ctemp = zeros(size(CNew));
            uniqueClass = unique(CNew);
            actualK = length(uniqueClass);
            for k = 1:actualK
            Ctemp(CNew==uniqueClass(k)) = k;
            end
            CNew =Ctemp;    
        end
    
        C = alignClusterings(Y(Y>1)-1,CNew);
        
        % Implement performance calculations
        confMat = confusionmat(Y(Y>1)-1,C);
    
    else
        % We consider the entire dataset
    
         C = alignClusterings(Y,C);
        
        % Implement performance calculations
        confMat = confusionmat(Y,C);
    
    end
    
    OA = sum(diag(confMat)/length(C)); 
    
    p=nansum(confMat,2)'*nansum(confMat)'/(nansum(nansum(confMat)))^2;
    kappa=(OA-p)/(1-p);
    
    ProdAcc = diag(confMat)./sum(confMat,2);
    mask = isnan(ProdAcc);
    ProdAcc(mask) = [];
    AA = mean(ProdAcc);
  
    
end