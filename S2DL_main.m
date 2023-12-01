%% Specify Dataset

clc
clear
 
profile off;
profile on;

prompt = 'Which dataset? \n 1) Indian Pines (Corrected) \n 2) Salinas (Corrected) \n 3) Salinas A (Corrected) \n';
DataSelected = input(prompt);
if DataSelected > 3 || DataSelected < 1
    disp('Incorrect prompt input. Please enter one of [1:3].')
end

datasets = {'IndianPines', 'Salinas', 'SalinasA'};

[X,M,N,D,HSI,GT,Y,~,~] = loadHSI(datasets{DataSelected});

[Idx_NN, Dist_NN] = knnsearch(X, X, 'K', 51);
Idx_NN = Idx_NN(:,2:51);
Dist_NN = Dist_NN(:,2:51);

% Set Default parameters
Hyperparameters.SpatialParams.ImageSize = [M,N];
Hyperparameters.NEigs = 10;
Hyperparameters.NumDtNeighbors = 200;
Hyperparameters.Beta = 2;
Hyperparameters.Tau = 10^(-5);
Hyperparameters.K_Known = length(unique(Y)); 
Hyperparameters.Tolerance = 1e-8;

clc
profile off
disp('Dataset Preloaded.')

%%
for nk = [1,3,4,5,6]
Hyperparameters.LocalBackbones = 1;
Hyperparameters.nk = nk;

for p = 100:100:1500
    
map = seg_ERS(HSI,0,p);
spSeg = double(map);
numSuperpixels = length(unique(spSeg));
Hyperparameters.Superpixel.map = spSeg;
Hyperparameters.Superpixel.num = numSuperpixels;
    
currentPerf = 0;
maxSum = NaN;

for l = 1:30
    NNs = 10:10:50;
    prctiles = 5:5:100;
    numReplicates = 1;

    OAs     = NaN*zeros(length(NNs), length(prctiles), numReplicates);
    kappas  = NaN*zeros(length(NNs), length(prctiles), numReplicates);
    AAs     = NaN*zeros(length(NNs), length(prctiles), numReplicates);
    Cs      = zeros(M*N,length(NNs), length(prctiles), numReplicates);

    disp('Set hyperparameter grid to begin grid search')
    Hyperparameters.SpatialParams.SpatialRadius = l;
    
    for i = 1:length(NNs)
        for j = 1:length(prctiles)
            for k = 1:numReplicates

                Hyperparameters.DiffusionNN = NNs(i);
                Hyperparameters.DensityNN = NNs(i); % must be ≤ 1000
                Hyperparameters.Sigma0 = prctile(Dist_NN(Dist_NN(:,1:NNs(i))>0), prctiles(j), 'all');

                density = KDE_large(Dist_NN, Hyperparameters);
                Clusterings = S2DL(X, density, Hyperparameters);

                if isfield(Clusterings,'Labels')
                    [ OAs(i,j,k), kappas(i,j,k), tIdx, ~, ~, AAs(i,j,k)] = calcPerformance(Y, Clusterings, ~strcmp('JasperRidge', datasets{DataSelected}));
                    C =  Clusterings.Labels(:,tIdx);
                    Cs(:,i,j,k) = C;
                end

                disp(['S2DL: '])
                disp([i/length(NNs), j/length(prctiles), k/numReplicates, currentPerf])

            end
            currentPerf = [max(nanmean(OAs,3),[],'all'), max(nanmean(kappas,3),[],'all'), max(nanmean(AAs,3),[],'all')];
        end 
        [n1,n2] = size(nanmean(OAs,3));
        [maxSum, k] = max(reshape(nanmean(OAs+kappas+AAs,3),n1*n2,1));
        [l,j] = ind2sub(size(mean(OAs+kappas+AAs,3)), k);
        save(strcat('S2DL', datasets{DataSelected}, num2str(numSuperpixels), 'SP', num2str(Hyperparameters.nk),'R', num2str(Hyperparameters.SpatialParams.SpatialRadius)),  'OAs', 'kappas', 'AAs', 'Cs', 'NNs', 'prctiles', 'numReplicates', 'maxSum')

    end
end
end
end
