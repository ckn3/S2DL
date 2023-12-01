function [V, VI_tot, t] = totalVI_minimization(Ct, Kt,  stability_weight)

[n, n_samp] = size(Ct);

J = find(and(Kt<n/2, Kt>1)); % Time samples during which a nontrivial clusering is extracted
notJ = setxor((1:n_samp), J);
n_J = length(J);
V = zeros(n_J); % To be the VI(Cs, Ct) matrix for nontrivial clusterings
VI_tot = zeros(n_samp,1);

if n_J > 0 % There is a nontrivial clustering of X.
    for i = 1:n_J
        for j = 1:n_J
            V(i,j) = VI(Ct(:,J(i)), Ct(:,J(j)));
        end
    end

    % Calculate Total VI
    if nargin == 2 % No stability weights. Each clustering is treated as the extracted clustering at time t.

        VI_tot(J) = sum(V);
        VI_tot(notJ) = NaN;

        [~,t] = min(VI_tot);

    elseif nargin ==3 % Stability weights

        VI_tot(J) = nansum(V,2).*stability_weight(J);  
        VI_tot(notJ) = NaN;

        [~,t] = min(VI_tot);

    end
    
else
    VI_tot(:) = NaN;
    t = 1;
end



