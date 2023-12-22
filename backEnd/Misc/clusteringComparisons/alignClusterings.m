function newLabels = alignClusterings(Y,Labels)
% Input:
%
% -Y: another assignment of labels of the same dataset as Y 
% -Labels: one assignment of labels of a dataset
%
% Output:
%
% -NewLabels: new version of Labels, designed to maximize correspondence with Y

% Allocate memory
n = length(Y);
newLabels=zeros(n,1);

% Find overlap between clusters   
K = length(unique(Y)); 

overlap=zeros(K,K); % ijth entry reflects the number of points in cluster i in 


for i=1:K
    for j=1:K
        overlap(i,j) = sum(and(Labels==i, Y==j));
    end
end

overlap=overlap';
 
% Run the Hungarian algorithm
cost=repmat(max(overlap,[],2),[1,K])-overlap;
idxOptimal = munkres(cost); % optimal indices 

% convert new Labels to 
for k=1:K
    newLabels(Labels==k) = find(idxOptimal(:,k));
end
end
