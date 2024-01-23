function Hyperparameters = loadHyperparameters(HSI, HSIName, AlgName)

if strcmp(HSIName, 'Synthetic Data') % Not an HSI 
    [n,D] = size(HSI);
    M=n; N=1;
    X = HSI;
else
    [M,N,D] = size(HSI);
    X = reshape(HSI,M*N,D);
end


if strcmp(AlgName, 'D-VIC')

    [~,Dist_NN] = knnsearch(X,X,'K', 1000);
    if strcmp(HSIName, 'Synthetic Data')
        NN  = 320;
        pct = 10.5;
        K = 3;
    elseif strcmp(HSIName, 'SalinasA')
        NN  = 30;
        pct = 96.2105263157895;
        K = 6;
    elseif strcmp(HSIName, 'Salinas')
        NN  = 900;
        pct = 20;
        K = 16;
    elseif strcmp(HSIName, 'JasperRidge')
        NN  = 20;
        pct = 92.758620689655170;
        K = 4;
    elseif strcmp(HSIName, 'IndianPines')
        NN = 80;
        pct = 100;
        K = 16;
    end

    % Set Default parameters
    Hyperparameters.SpatialParams.ImageSize = [M,N];
    Hyperparameters.NEigs = 10;
    Hyperparameters.NumDtNeighbors = 200;
    Hyperparameters.Beta = 2;
    Hyperparameters.Tau = 10^(-5);
    Hyperparameters.Tolerance = 1e-8;
    Hyperparameters.K_Known = K; 
    Hyperparameters.EndmemberParams.Algorithm = 'ManyAVMAX';
    Hyperparameters.EndmemberParams.NumReplicates = 100;
    Hyperparameters.EndmemberParams.K = hysime(X');
    Hyperparameters.DiffusionNN = NN;
    Hyperparameters.DensityNN = NN; % must be ≤ 1000
    Hyperparameters.Sigma0 = prctile(Dist_NN(Dist_NN>0), pct, 'all');

elseif strcmp(AlgName, 'S2DL')

    if strcmp(HSIName, 'Synthetic Data')
        R = 0;
        p = 0;
        nk = 0;
        NN  = 0;
        pct = 0;
        K = 3;
    elseif strcmp(HSIName, 'SalinasA')
        R = 7;
        p = 1000;
        nk = 5;
        NN  = 50;
        pct = 10;
        K = 6;
    elseif strcmp(HSIName, 'Salinas')
        R = 11;
        p = 1500;
        nk = 3;
        NN  = 10;
        pct = 70;
        K = 16;
    elseif strcmp(HSIName, 'JasperRidge')
        R = 0;
        p = 0;
        nk = 0;
        NN  = 0;
        pct = 0;
        K = 4;
    elseif strcmp(HSIName, 'IndianPines')
        R = 9;
        p = 1000;
        nk = 5;
        NN = 20;
        pct = 55;
        K = 16;
    end

    % Set Default parameters
    [~,Dist_NN] = knnsearch(X,X,'K', NN+1);
    Dist_NN(:,1) = [];
    Hyperparameters.LocalBackbones = 1;
    Hyperparameters.nk = nk;
    Hyperparameters.p = p;
    Hyperparameters.SpatialParams.SpatialRadius = R;
    Hyperparameters.SpatialParams.ImageSize = [M,N];
    Hyperparameters.NEigs = 10;
    Hyperparameters.NumDtNeighbors = 200;
    Hyperparameters.Beta = 2;
    Hyperparameters.Tau = 10^(-5);
    Hyperparameters.Tolerance = 1e-8;
    Hyperparameters.K_Known = K; 
    Hyperparameters.DiffusionNN = NN;
    Hyperparameters.DensityNN = NN; % must be ≤ 1000
    Hyperparameters.Sigma0 = prctile(Dist_NN(Dist_NN>0), pct, 'all');

elseif strcmp(AlgName, 'LUND')

    [~,Dist_NN] = knnsearch(X,X,'K', 1000);
    if strcmp(HSIName, 'Synthetic Data')
        NN  = 140;
        pct = 0.5;
        K = 3;
    elseif strcmp(HSIName, 'SalinasA')
        NN  = 40;
        pct = 5;
        K = 6;
    elseif strcmp(HSIName, 'Salinas')
        NN  = 900;
        pct = 85;
        K = 16;
    elseif strcmp(HSIName, 'JasperRidge')
        NN  = 40;
        pct = 75;
        K = 4;
    elseif strcmp(HSIName, 'IndianPines')
        NN = 40;
        pct = 75;
        K = 16;
    end

    % Set Default parameters
    Hyperparameters.SpatialParams.ImageSize = [M,N];
    Hyperparameters.NEigs = 10;
    Hyperparameters.NumDtNeighbors = 200;
    Hyperparameters.Beta = 2;
    Hyperparameters.Tau = 10^(-5);
    Hyperparameters.Tolerance = 1e-8;
    Hyperparameters.K_Known = K; 
    Hyperparameters.DiffusionNN = NN;
    Hyperparameters.DensityNN = NN; % must be ≤ 1000
    Hyperparameters.Sigma0 = prctile(Dist_NN(Dist_NN>0), pct, 'all');

elseif strcmp(AlgName, 'SRDL')

    [~,Dist_NN] = knnsearch(X,X,'K', 50);
    if strcmp(HSIName, 'SalinasA')
        NN  = 0;
        pct = 0;
        R = 0;
        K = 6;
        
    elseif strcmp(HSIName, 'Salinas')
        NN  = 10;
        pct = 5;
        R = 5;
        K = 16;
        
    elseif strcmp(HSIName, 'JasperRidge')
        NN  = 0;
        pct = 0;
        R = 0;
        K = 4;
        
    elseif strcmp(HSIName, 'IndianPines')
        NN = 10;
        pct = 15;
        R = 4;
        K = 16;
        Hyperparameters.SpatialParams.ConsensusSpatialRadius = 2;

    elseif strcmp(HSIName, 'WHU')
        NN = 0;
        pct = 0;
        R = 0;
        K = 9;
        Hyperparameters.SpatialParams.ConsensusSpatialRadius = 2;
    end
    
    % Set Default parameters
    Hyperparameters.SpatialParams.ImageSize = [M,N];
    Hyperparameters.NEigs = 10;
    Hyperparameters.NumDtNeighbors = 200;
    Hyperparameters.Beta = 2;
    Hyperparameters.Tau = 10^(-5);
    Hyperparameters.Tolerance = 1e-8;
    Hyperparameters.SpatialParams.GraphSpatialRadius = R;
    Hyperparameters.SpatialParams.ConsensusSpatialRadius = R;
    Hyperparameters.K_Known = K; 
    Hyperparameters.DiffusionNN = NN;
    Hyperparameters.DensityNN = NN; % must be ≤ 1000
    Hyperparameters.Sigma0 = prctile(Dist_NN(Dist_NN>0), pct, 'all');

elseif strcmp(AlgName, 'DSIRC')

    [~,Dist_NN] = knnsearch(X,X,'K', 1000);
    if strcmp(HSIName, 'Synthetic Data')
        NN  = 140;
        pct = 0.5;
        K = 3;
    elseif strcmp(HSIName, 'SalinasA')
        NN  = 60;
        pct = 92.8947368421053;
        K = 6;
    elseif strcmp(HSIName, 'Salinas')
        NN  = 30;
        pct = 100;
        K = 16;
    elseif strcmp(HSIName, 'JasperRidge')
        NN  = 0;
        pct = 0;
        K = 4;
    elseif strcmp(HSIName, 'IndianPines')
        NN = 30;
        pct = 72.0526315789474;
        K = 16;
    end

    % Set Default parameters
    Hyperparameters.SpatialParams.ImageSize = [M,N];
    Hyperparameters.NEigs = 10;
    Hyperparameters.NumDtNeighbors = 200;
    Hyperparameters.Beta = 2;
    Hyperparameters.Tau = 10^(-5);
    Hyperparameters.Tolerance = 1e-8;
    Hyperparameters.K_Known = K; 
    Hyperparameters.EndmemberParams.Algorithm = 'ManyAVMAX';
    Hyperparameters.EndmemberParams.NumReplicates = 100;
    Hyperparameters.EndmemberParams.K = hysime(X');
    Hyperparameters.DiffusionNN = NN;
    Hyperparameters.DensityNN = NN; % must be ≤ 1000
    Hyperparameters.Sigma0 = prctile(Dist_NN(Dist_NN>0), pct, 'all');
    
elseif strcmp(AlgName, 'SC')

    if strcmp(HSIName, 'Synthetic Data')
        error('This algorithm is not supported for Synthetic Dataset')
    elseif strcmp(HSIName, 'SalinasA')
        NN  = 20;
        K = 6;
    elseif strcmp(HSIName, 'Salinas')
        NN  = 500;
        K = 16;
    elseif strcmp(HSIName, 'JasperRidge')
        NN  = 10;
        K = 4;
    elseif strcmp(HSIName, 'IndianPines')
        NN = 100;
        K = 16;
    end

    Hyperparameters.DiffusionNN = NN;
    Hyperparameters.NEigs = min(K+1, 10); 
    Hyperparameters.SpatialParams.ImageSize = [M,N];    

elseif strcmp(AlgName, 'PGDPC')

    if strcmp(HSIName, 'Synthetic Data')
        NN  = 0;
        K = 3;
    elseif strcmp(HSIName, 'SalinasA')
        NN  = 30;
        K = 6;
    elseif strcmp(HSIName, 'Salinas')
        NN  = 0;
        K = 16;
    elseif strcmp(HSIName, 'JasperRidge')
        NN  = 0;
        K = 4;
    elseif strcmp(HSIName, 'IndianPines')
        NN = 100;
        K = 16;
    end

    % Set Default parameters
    Hyperparameters.K_Known = K; 
    Hyperparameters.NN = NN;
    
elseif strcmp(AlgName, 'SPGDPC')

    if strcmp(HSIName, 'Synthetic Data')
        NN  = 0;
        supN = 0;
        sigma = 0;
        K = 3;
    elseif strcmp(HSIName, 'SalinasA')
        NN  = 10;
        supN = 1200;
        sigma = 1.25892541179417;
        K = 6;
    elseif strcmp(HSIName, 'Salinas')
        NN  = 20;
        supN = 800;
        sigma = 10;
        K = 16;
    elseif strcmp(HSIName, 'JasperRidge')
        NN  = 0;
        supN = 0;
        sigma = 0;
        K = 4;
    elseif strcmp(HSIName, 'IndianPines')
        NN = 40;
        supN = 1100;
        sigma = 15.8489319246111;
        K = 16;
    end

    % Set Default parameters
    Hyperparameters.K_Known = K; 
    Hyperparameters.NN = NN;
    Hyperparameters.supN = supN;
    Hyperparameters.sigma = sigma;
    
elseif strcmp(AlgName, 'SCI')

    if strcmp(HSIName, 'Synthetic Data')
        beta  = 0;
        pct = 0;
        K = 3;
    elseif strcmp(HSIName, 'SalinasA')
        beta  = 0.01;
        pct = 15;
        K = 6;
    elseif strcmp(HSIName, 'Salinas')
        beta  = 0;
        pct = 0;
        K = 16;
    elseif strcmp(HSIName, 'JasperRidge')
        beta  = 0;
        pct = 0;
        K = 4;
    elseif strcmp(HSIName, 'IndianPines')
        beta = 0.630957344480193;
        pct = 5;
        K = 16;
    end

    % Set Default parameters
    Hyperparameters.K_Known = K; 
    Hyperparameters.beta = beta;
    Hyperparameters.prctiles = pct;

elseif strcmp(AlgName, 'DLSS')

    [~,Dist_NN] = knnsearch(X,X,'K', 900);
    if strcmp(HSIName, 'SalinasA')
        NN  = 0;
        pct = 0;
        R = 0;
        K = 6;
        
    elseif strcmp(HSIName, 'Salinas')
        NN  = 0;
        pct = 0;
        R = 0;
        K = 16;
        
    elseif strcmp(HSIName, 'JasperRidge')
        NN  = 0;
        pct = 0;
        R = 0;
        K = 4;
        
    elseif strcmp(HSIName, 'IndianPines')
        NN = 0;
        pct = 0;
        R = 0;
        K = 16;

    elseif strcmp(HSIName, 'WHU')
        NN = 0;
        pct = 0;
        R = 0;
        K = 9;
        
    end
    
    % Set Default parameters
    Hyperparameters.SpatialParams.ImageSize = [M,N];
    Hyperparameters.NEigs = 10;
    Hyperparameters.NumDtNeighbors = 200;
    Hyperparameters.Beta = 2;
    Hyperparameters.Tau = 10^(-5);
    Hyperparameters.Tolerance = 1e-8;
    Hyperparameters.SpatialParams.ConsensusSpatialRadius = R;
    Hyperparameters.K_Known = K; 
    Hyperparameters.DiffusionNN = NN;
    Hyperparameters.DensityNN = NN; % must be ≤ 1000
    Hyperparameters.Sigma0 = prctile(Dist_NN(Dist_NN>0), pct, 'all');
    
else 
    Hyperparameters = [];
end
