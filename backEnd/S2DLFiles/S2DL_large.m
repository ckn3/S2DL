function [C, K, Dt] = S2DL_large(X, Hyperparameters, t, G, p)
%{
 - This function produces a structure with multiscale clusterings produced
   with the LUND algorithm, presented in the following paper. 
        - Maggioni, M., J.M. Murphy. Learning by Unsupervised Nonlinear 
          Diffusion. Journal of Machine Learning Research, 20(160), 
          pp. 1-56. 2019.
    
   and analyzed further in the following papers:
    
        - Murphy, J. M., & Polk, S. L. (2022). A multiscale environment for 
          learning by diffusion. Applied and Computational Harmonic 
          Analysis, 57, 58-100.
        - Polk, S. L., & Murphy, J. M. (2021, July). Multiscale Clustering 
          of Hyperspectral Images Through Spectral-Spatial Diffusion 
          Geometry. In 2021 IEEE International Geoscience and Remote 
          Sensing Symposium IGARSS (pp. 4688-4691). IEEE.
        - Polk, S. L., Cui, K., Plemmons, R. J., and Murphy, J. M., (2022). 
          Diffusion and Volume Maximization-Based Clustering of Highly 
          Mixed Hyperspectral Images. (In Review).
Inputs: X:                      Data matrix.
        Hyperparameters:        Optional structure with graph parameters
                                with required fields:  
        t:                      Diffusion time parameter.
        G:                      Graph structure computed using  
                                'extract_graph_large.m' 
        p:                      Kernel Density Estimator.
        time_idx:                      Kernel Density Estimator.
Output: 
            - C:                n x 1 vector storing the LUND clustering 
                                of X at time t.
            - K:                Scalar, number of clusters in C.
                                clusters in the Labels(:,t) clustering.
            - Dt:               n x 1 matrix storing \mathcal{D}_t(x). 
© 2022 Sam L Polk, Tufts University. 
email: samuel.polk@tufts.edu
%}  

if ~isfield(Hyperparameters, 'NEigs')
    Hyperparameters.NEigs = size(G.EigenVecs,2);
end
if ~isfield(Hyperparameters, '')
    Hyperparameters.DtNNs = 100;
end

[n,~] = size(X);

if Hyperparameters.LocalBackbones == 2
    % Calculate diffusion map using the spatial agnostic graph
    DiffusionMap = zeros(n,Hyperparameters.NEigs);
    for l = 1:size(DiffusionMap,2)
        DiffusionMap(:,l) = G.G1.EigenVecs(:,l).*(G.G1.EigenVals(l).^t);
    end
else
    % Calculate diffusion map
    DiffusionMap = zeros(n,Hyperparameters.NEigs);
    for l = 1:size(DiffusionMap,2)
        DiffusionMap(:,l) = G.EigenVecs(:,l).*(G.EigenVals(l).^t);
    end
end
 
% Compute Hyperparameters.NumDtNeighbors Dt-nearest neighbors.
[IdxNN, D] = knnsearch(DiffusionMap, DiffusionMap, 'K', Hyperparameters.NumDtNeighbors);

% compute rho_t(x), stored as rt
rt = zeros(n,1);
% p = p./sum(p);
p_max = max(p);
for i=1:n
   
    if p(i) == p_max
        rt(i) = max(pdist2(DiffusionMap(i,:),DiffusionMap));
    else
        
        idces =  find(p(IdxNN(i,:))>p(i));
        
        if ~isempty(idces)
            % In this case, at least one of the Hyperparameters.NumDtNeighbors Dt-nearest neighbors
            % of X(i,:) is also higher density, so we have already calculated rho_t(X(i,:)).
            rt(i) = D(i,idces(1));
            
        else
            % In this case, none of the first Hyperparameters.NumDtNeighbors Dt-nearest neighbors of
            % X(i,:) are also higher density. So, we do the full search.
            rt(i) = min(pdist2(DiffusionMap(i,:),DiffusionMap(p>p(i),:)));
        end
    end
end
 
% Extract Dt(x) and sort in descending order
Dt = rt.*p;
[~, m_sorting] = sort(Dt,'descend');

% Determine K based on the ratio of sorted Dt(x_{m_k}). 
if isfield(Hyperparameters, 'K_Known')
    K = Hyperparameters.K_Known;
else
    [~, K] = max(Dt(m_sorting(2:ceil(n/2)-1))./Dt(m_sorting(3:ceil(n/2))));
    K=K+1;
end

if K == 1
    C = ones(n,1);
else
    
    idx = 1:n;
    C = zeros(n,1);
    % Label modes
    % C(m_sorting(1:K)) = 1:K;
    
    % Make sure each superpixel can have at most 1 mode.
    % m1 indicates the superpixel that it belongs to.
    m1=floor((m_sorting-1)/Hyperparameters.nk)+1;
    if length(unique(m1(1:K)))<K
        % loop until the first few indices of m1 contain K unique values.
        i=0;
        while length(unique(m1(1:K+i)))<K
            i=i+1;
        end
        m1=m1(1:K+i);
        C(m_sorting(1)) = 1;
        lb=2;
        for j=2:length(m1)
            if length(unique(m1(1:j-1)))<length(unique(m1(1:j)))
                C(m_sorting(j)) = lb;
                lb=lb+1;
            else
                % I = find(m1==m1(j),1);
                % C(m_sorting(j)) = C(m_sorting(I));
                C(m_sorting(j)) = 0;
            end
        end
    else
        C(m_sorting(1:K)) = 1:K;
    end
    
    % Forming local backbones, i.e., giving nearest neighbors of modes the
    % same label as these modes.
    % DistLB = D(m_sorting(1:K),2:Hyperparameters.DiffusionNN);
    if Hyperparameters.LocalBackbones == 1
        % IdxLB contains the spatial-regularized NNs of modes.
        IdxLB = IdxNN(m_sorting(C(m_sorting)>0),2:Hyperparameters.DiffusionNN); % Hyperparameters.DiffusionNN max(10,round(0.1*n/K))
        for j = 1:Hyperparameters.DiffusionNN-1
            % if the column are unlabeled, label them from 1 to K.
            if sum(C(IdxLB(:,j))) == 0
                C(IdxLB(:,j)) = 1:K;
            % otherwise, only label those unlabeled NNs of modes.
            else
                for k = 1:K
                    if C(IdxLB(k,j)) == 0
                        C(IdxLB(k,j)) = k;
                    end
                end
            end
        end
%     elseif Hyperparameters.LocalBackbones == 2
%         % Construct another Diffusion map for the spatial agnostic graph.
%         DiffusionMap1 = zeros(n,Hyperparameters.NEigs);
%         for l = 1:size(DiffusionMap1,2)
%             DiffusionMap1(:,l) = G.G1.EigenVecs(:,l).*(G.G1.EigenVals(l).^t);
%         end
%         
%         % Search the NNs in the spatial agnostic graph and use local backbone.
%         [IdxNN1, ~] = knnsearch(DiffusionMap1, DiffusionMap1, 'K', Hyperparameters.NumDtNeighbors);
%         IdxLB = IdxNN1(m_sorting(C(m_sorting)>0),2:Hyperparameters.DiffusionNN);
%         for j = 1:Hyperparameters.DiffusionNN-1
%             % if the column are unlabeled, label them from 1 to K.
%             if sum(C(IdxLB(:,j))) == 0
%                 C(IdxLB(:,j)) = 1:K;
%             % otherwise, only label those unlabeled NNs of modes.
%             else
%                 for k = 1:K
%                     if C(IdxLB(k,j)) == 0
%                         C(IdxLB(k,j)) = k;
%                     end
%                 end
%             end
%         end
    end

    % Label non-modal points according to the label of their Dt-nearest
    % neighbor of higher density that is already labeled.
    [~,l_sorting] = sort(p,'descend');
     
    for j = 1:n
        i = l_sorting(j);
        if C(i)==0 % unlabeled point
            
            NNs = IdxNN(i,:);
            idces = find(and(C(NNs)>0, p(NNs)>p(i))); % Labeled, higher-density points in NNs
            
            if isempty(idces)
                % None of the Dt-nearest neighbors are also higher-density
                % & labeled. So, we do a full search.
                candidates = idx(and(C>0, p>=p(i))); % All labeled points of higher density.
                if isempty(candidates)
                    disp([])
                end
                
                [~,temp_idx] = min(pdist2(DiffusionMap(i,:), DiffusionMap(candidates,:)));
                C(i) = C(candidates(temp_idx));
            else
                % At least one of the Dt-nearest neighbors is higher
                % density & labeled. So, we pick the closest point.
                C(i) = C(NNs(idces(1)));                
            end
        end
    end
    
end 