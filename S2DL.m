function [Clusterings, runtimes] = S2DL(X, density, Hyperparameters)

purityMode = 0;
k = Hyperparameters.nk;

tic;
if purityMode == 0
    % Use density
    zeta = density;
    spSeg = Hyperparameters.Superpixel.map;
    numSuperpixels = Hyperparameters.Superpixel.num;
    spFlat = spSeg(:);
    idx = zeros(numSuperpixels*k,1);

    for i = 1:numSuperpixels

        % create mask for the i-1st superpixel
        mask_z = zeros(length(spFlat),1);
        mask_z(spSeg==(i-1)) = zeta(spSeg==(i-1));

        % save the corresponding index and zeta
        [~,order] = sort(mask_z,'descend');
        idx((k*(i-1)+1):k*i) = order(1:k);
    end

    Xsp = X(idx,:);
    zetaSp = zeta(idx);
    
elseif purityMode == 1
    kk = 2*k-1;
    % Use density and consider distance
    zeta = density;
    spSeg = Hyperparameters.Superpixel.map;
    numSuperpixels = Hyperparameters.Superpixel.num;
    spFlat = spSeg(:);
    idx = zeros(numSuperpixels*k,1);

    for i = 1:numSuperpixels

        % create mask for the i-1st superpixel
        mask_z = zeros(length(spFlat),1);
        mask_z(spSeg==(i-1)) = zeta(spSeg==(i-1));

        % save the corresponding index and zeta
        [~,order] = sort(mask_z,'descend');
        order = order(1:kk);
        % calculate the coords of these pixels
        xcd = mod(order-1,kk)+1;
        ycd = floor((order-1)/5)+1;
        % find the pixels with largest (1/dist * density) to the density
        % maximizer of each superpixel
        [~,order1] = sort(sqrt(pdist2([xcd(1),ycd(1)],[xcd(2:kk),ycd(2:kk)])').*zeta(order(2:kk)),'descend');
        idx((k*(i-1)+1):k*i) = [order(1); order(order1(1:(k-1))+1)];
    end

    Xsp = X(idx,:);
    zetaSp = zeta(idx);
    
elseif purityMode == 2
    % Use density and purity
    [pixelPurity, ~, ~] = compute_purity(X,Hyperparameters);
    zeta = harmmean([density./max(density), pixelPurity./max(pixelPurity)],2);

    spSeg = Hyperparameters.Superpixel.map;
    numSuperpixels = Hyperparameters.Superpixel.num;
    spFlat = spSeg(:);
    idx = zeros(numSuperpixels*k,1);

    for i = 1:numSuperpixels

        % create mask for the i-1st superpixel
        mask_z = zeros(length(spFlat),1);
        mask_z(spSeg==(i-1)) = zeta(spSeg==(i-1));

        % save the corresponding index and zeta
        [~,order] = sort(mask_z,'descend');
        idx((k*(i-1)+1):k*i) = order(1:k);
    end

    Xsp = X(idx,:);
    zetaSp = zeta(idx);
    
else
    % Use density then purity
    spSeg = Hyperparameters.Superpixel.map;
    numSuperpixels = Hyperparameters.Superpixel.num;
    spFlat = spSeg(:);
    idx = zeros(numSuperpixels*k,1);

    for i = 1:numSuperpixels

        % create mask for the i-1st superpixel
        mask_z = zeros(length(spFlat),1);
        mask_z(spSeg==(i-1)) = density(spSeg==(i-1));

        % save the corresponding index and zeta
        [~,order] = sort(mask_z,'descend');
        idx((k*(i-1)+1):k*i) = order(1:k);
    end

    Xsp = X(idx,:);
    
    [pixelPurity, ~, ~] = compute_purity(Xsp,Hyperparameters);
    zetaSp = harmmean([density(idx)./max(density(idx)), pixelPurity./max(pixelPurity)],2);
end

Hyperparameters.Superpixel.idx=idx;

[G,~] = extract_graph_superpixel(Xsp, Hyperparameters);
runtimeTemp1 = toc;

if G.EigenVals(2)<1
    [Clusterings, runtimes] = MS2DL_large(Xsp, Hyperparameters, G, zetaSp);
else
    Clusterings = NaN;
    runtimes = 0;
end  

tic;
[M,N] = size(spSeg);
pred = zeros(M,N);
if isfield(Clusterings,'K')
    for t = 1:length(Clusterings.K)
        temp = Clusterings.Labels(:,t);
        for ii = 1:numSuperpixels
            mask = (spSeg==(ii-1));
            pred(mask)=mode(temp((k*(ii-1)+1):(k*ii)));
        end
        Label(:,t) = reshape(pred,M*N,1);
    end
    Clusterings.Labels = Label;
end
runtimeTemp2 = toc;

runtimes = runtimes + runtimeTemp1+runtimeTemp2;

end