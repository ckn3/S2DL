function [NMI, MI, Hu, Hv] = nmi(U,V)
%{
Calculates the normalized mutual information between, mutual information 
between, and entropy of clusterings U and V of a dataset X. 

Inputs:     U:  nx1 vector encoding first clustering of X.
            V:  nx1 vector encoding second clustering of X.

Outputs:    NMI:    The normalized mutual information between U and V. 
            MI:     The mutual information between U and V. 
            Hu:     The entropy of U.
            Hv:     The entropy of V.
%}

assert(numel(U) == numel(V)); %  Breaks script if not right dimensions.

n = numel(U); % size of dataset

% converts to the right size
U = reshape(U,1,n);
V = reshape(V,1,n); 

% Makes sure that U and V start at k=1. 
l = min(min(U),min(V));
U = U-l+1;
V = V-l+1;
k_max = max(max(U),max(V));

%Get distributions
idx = 1:n;
Mu = sparse(idx,U,1,n,k_max,n);
Mv = sparse(idx,V,1,n,k_max,n); 
Puv = nonzeros(Mu'*Mv/n);   % Joint distribution of U and V.
Huv = -dot(Puv,log2(Puv));  % Joint Entropy of U and V

Pu = nonzeros(mean(Mu,1));  % Distribution of U
Pv = nonzeros(mean(Mv,1));  % Distribution of V

% Calculate Entropy of U and V:
Hu = -dot(Pu,log2(Pu));
Hv = -dot(Pv,log2(Pv));

% Calculate MI between U and V
MI = Hu + Hv - Huv;

% Calculate NMI(U,V)
NMI = sqrt((MI/Hu)*(MI/Hv));
NMI = MI/(Hu+Hv)*2;
NMI = max(0,NMI);
