function W = spatial_weight_matrix(X,Hyperparameters)
%{
 - This function computes the weight matrix for a modified graph that 
   directly incorporates spatial geometry into graph diffusion. See paper 
   below for more.

        - Polk, Sam L. and Murphy James M., 2021. Multiscale Spectral-
          Spatial Diffusion Geometry for Hyperspectral Image Clustering. 
          (In Review)

Inputs: X:                      (M*N)xD Data matrix .
        Hyperparameters:        Structure with graph parameters with the 
                                required fields: 
            - SpatialParams:    Stores the dimensions of the original image
            - DiffusionNN:      Number of nearest neighbors in KNN graph.
            - WeightType:       Equal to either 'adjesency' or 'gaussian' (Optional).
            - Sigma:            If WeightType == 'gaussian', then diffusion scale parameter Sigma>0 required.

Output: W:                  Weight matrix for graph.

© 2021 Sam L Polk, Tufts University. 
email: samuel.polk@tufts.edu
%}


% If spatial information is included in hyperparameters structure, we
% incorporate that into the diffusion process.
R = Hyperparameters.SpatialParams.SpatialRadius;
M = Hyperparameters.SpatialParams.ImageSize(1); 
N = Hyperparameters.SpatialParams.ImageSize(2);
n = size(X,1); % Number of pixels in X
NN = Hyperparameters.DiffusionNN;

if strcmp(Hyperparameters.WeightType, 'gaussian')
    sigma = Hyperparameters.Sigma;
end 

idx_row = zeros(n*NN,1);  % Rows where nonzero elements appear in W
idx_col = zeros(n*NN,1);  % Columns where nonzero elements appear in W
vals    = zeros(n*NN,1);  % Variable to store edge weights

ct = 1;
% neighbors = sur(Hyperparameters.Superpixel.map);
idx_sp = Hyperparameters.Superpixel.idx;
temp=1:length(idx_sp);
[i_t,j_t]=ind2sub([M,N],idx_sp);

for idx = 1:n
    
    % Find the spatial nearest neighbors of the point X(idx,:)
    [i,j] = ind2sub([M,N],idx_sp(idx));
    % Find indices of pixel neighbors of pixel idx(j)
%     temp=[];
%     nbh = unique(neighbors(1,neighbors(2,:)==floor((idx-1)/Hyperparameters.nk)));
%     for k=1:Hyperparameters.nk
%         temp = [temp, (nbh+1)*Hyperparameters.nk-k+1];
%     end
%     temp=sort(temp);
%     [i_t,j_t]=ind2sub([M,N],idx_sp(temp));
    flag = abs(i_t-i)<=R & abs(j_t-j)<=R;
    NN_Idx=temp(flag);
    % [i,j] = ind2sub([M,N], idx);
    % NN_Idx=FindNeighbors([i,j], R, M, N); % indices of the spatial nearest neighbors of ij pixel of X.
    NN_Count=length(NN_Idx); % Number of spatial nearest neighbors of X(idx,:);
    
    Xi_rep = repmat(X(idx,:),NN_Count,1); % Spectra of ij pixel of X, repeated in each row.
    X_NNs = X(NN_Idx,:); % Spectra of spatial nearest neighbors of ij pixel
        
    DistTemp=sqrt(sum((Xi_rep - X_NNs).^2,2)); % distance between X(idx,:) and its spatial nearest neighbors
    [NN_Dists,Idx]=sort(DistTemp,'ascend');
    NN_Idx=NN_Idx(Idx); 
    
    nEdge = min(NN_Count, NN); 
    
    idx_row(ct:ct+nEdge-1) = idx;
    idx_col(ct:ct+nEdge-1) = NN_Idx(1:nEdge);
    
    if strcmp(Hyperparameters.WeightType, 'adjesency')
        vals(ct:ct+nEdge-1) = 1; 
    else
        vals(ct:ct+nEdge-1) = exp(-(NN_Dists(1:nEdge).^2)./(sigma^2));
    end  

    ct = ct+nEdge;
    % disp(idx/n)
end

% Truncate values we didn't use.
idx_row(ct:end) = [];
idx_col(ct:end) = [];
vals(ct:end) = [];

W = sparse(idx_row, idx_col, vals);

if size(W,1)>size(W,2)
    W = 0.5*([W, zeros(n, n-size(W,2))] + [W, zeros(n, n-size(W,2))]');
else    
    W = adjacency((graph((W+W')./2)));   % Convert W to adjesency matrix. (W+W)./2 forces symmetry.
end

end

