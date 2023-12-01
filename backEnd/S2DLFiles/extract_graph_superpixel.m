function [G,W] = extract_graph_superpixel(X, Hyperparameters, Idx_NN, Dist_NN)

[n,~] = size(X); 
NN = Hyperparameters.DiffusionNN;
n_eigs = Hyperparameters.NEigs;


if ~isfield(Hyperparameters, 'WeightType')
    Hyperparameters.WeightType = 'adjesency';
end

if isfield(Hyperparameters.SpatialParams, 'SpatialRadius') 

    if strcmp(Hyperparameters.WeightType, 'adjesency')
        W = spatial_weight_matrix_sp(X,Hyperparameters);
    else
        if ~isfield(Hyperparameters, 'Sigma')
            Hyperparameters.Sigma = prctile(Dist_NN(:,1:NN),50,'all');
        end
        W = spatial_weight_matrix(X,Hyperparameters);
    end
    
    % Construct a spatial agnostic graph for local backbone
    if Hyperparameters.LocalBackbones == 2
        [Idx_NN, ~] = knnsearch(X,X,'K',51);
        ind_row = repmat((1:n)', 1,NN);  % row indices for nearest neighbors.
        ind_col = Idx_NN(:,1:NN);         % column indices for nearest neighbors.
        W1 = sparse(ind_row, ind_col, ones(n,NN)); % construct W
        if size(W1,1)>size(W1,2)
            W1 = 0.5*([W1, zeros(n, n-size(W1,2))] + [W1, zeros(n, n-size(W1,2))]');
        else    
            W1 = adjacency((graph((W1+W1')./2)));   % Convert W to adjesency matrix. (W+W)./2 forces symmetry.
        end
        W1 = (W1+W1')./2;
        [V1,D1, ~] = eigs(spdiags(1./sum(W1)',0,n,n)*W1, n_eigs, 'largestabs'); 
        [lambda1,idx1] = sort(diag(abs(D1)),'descend');
        lambda1(1) = 1;
        V1(:,1) = 1;
        V1 = real(V1(:,idx1));
        G1.EigenVecs = V1;
        G1.EigenVals = lambda1;
    end
    
else
    
    % Preallocate memory for sparse matrix calculation
    ind_row = repmat((1:n)', 1,NN);  % row indices for nearest neighbors.
    ind_col = Idx_NN(:,1:NN);         % column indices for nearest neighbors.
    if strcmp(Hyperparameters.WeightType, 'adjesency')

        W = sparse(ind_row, ind_col, ones(n,NN)); % construct W
        if size(W,1)>size(W,2)
            W = 0.5*([W, zeros(n, n-size(W,2))] + [W, zeros(n, n-size(W,2))]');
        else    
            W = adjacency((graph((W+W')./2)));   % Convert W to adjesency matrix. (W+W)./2 forces symmetry.
        end
        
    elseif strcmp(Hyperparameters.WeightType, 'gaussian')

        if ~isfield(Hyperparameters, 'Sigma')
            sigma = prctile(Dist_NN(:,1:NN),50,'all');
        else
            sigma = Hyperparameters.Sigma;
        end

        W = sparse(ind_row, ind_col, exp(-(Dist_NN(:,1:NN).^2)./(sigma^2))); % construct W
    end
end

W = (W+W')./2;   %  (W+W)./2 forces symmetry

% First n_eigs eigenpairs of transition matrix D^{-1}W:
[V,D, flag] = eigs(spdiags(1./sum(W)',0,n,n)*W, n_eigs, 'largestabs'); 
    
if flag
    % Didn't converge. Try again with larger subspace dimension.
    [V,D, flag] = eigs(spdiags(1./sum(W)',0,n,n)*W, n_eigs, 'largestabs',  'SubspaceDimension', max(4*n_eigs,40));  
end


if flag
    disp('Convergence Failed.')
    G = NaN;
else
    [lambda,idx] = sort(diag(abs(D)),'descend');
    lambda(1) = 1;
    V(:,1) = 1;
    V = real(V(:,idx));
    

    G.Hyperparameters = Hyperparameters;
    G.EigenVecs = V;
    G.StationaryDist = sum(W)./sum(W, 'all');
    G.EigenVals = lambda;
    
    % Save
    if Hyperparameters.LocalBackbones == 2
        G.G1 = G1;
    end
    
end


