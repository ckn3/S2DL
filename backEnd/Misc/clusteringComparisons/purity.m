function [ p ] = purity( pre, label )
%CALACCURANY Summary of this function goes here
%   Detailed explanation goes here
tab = crosstab(pre, label);
[Ral,Row] = max(tab,[],1);
[Cal,Col] = max(tab,[],2);
% for i = 1:size(Row,1)
for i = 1:size(Col,1)
    summ(i) = itersum( 0, tab, i, Col(i) );
end
for j = 1:size(Row,2)
    summ(i+j) = itersum( 0, tab, Row(j), j );
end
p = max(summ)./sum(sum(tab));
end

function [ msum ] = itersum( sum, table, row, col )
    if(size(table,2) == 1 || size(table,1) == 1 )
        sum = sum + table(row, col);
    else
        sum = sum + table(row, col);
        table(row,:) = [];
        table(:,col) = [];
        [mv,mr] = max(table,[],1);
        [~,mc] = max(mv,[],2);
        sum = itersum(sum, table, mr(mc), mc);
    end
    msum = sum;
end