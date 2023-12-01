function [vi, NMI] = VI(U,V)
%{
Calculates the variation of information and normalized mutual information
between clusterings U and V of a dataset X. 

Inputs:     U:  nx1 vector encoding first clustering of X.
            V:  nx1 vector encoding second clustering of X.

Outputs:    V:  The variation of information between U and V. 
            U:  The normalized mutual information between U and V. 

%}

assert(numel(U) == numel(V)); %  Breaks script if not right dimensions.
[NMI, MI, Hu, Hv] = nmi(U,V);

% calculate VI(U,V)
vi = Hu + Hv - 2*MI; 
