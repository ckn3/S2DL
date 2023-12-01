function p = KDE_large(Dist_NN, Hyperparameters)
%{
Inputs:     Dist_NN:            nxM matrix where Dist_NN(i,:) encodes the M 
                                nearest neighbors of X(i,:), sorted in 
                                ascending order.
            Hyperparameters:    Structure with 

Outputs:    p:                  Kernel density estimator. 
%}

% Extract hyperparameters
NN = Hyperparameters.DensityNN;
sigma0 = Hyperparameters.Sigma0;

% Calculate density
D_trunc = Dist_NN(:,1:NN);
p = sum(exp(-(D_trunc.^2)./(sigma0^2)),2);
p = p./sum(p);


% avoids p=0 by setting 0-value entries as the minimum value 
p(p == 0) = min(p(p>0)); 