function [pvalues,stats_values] = compute_pvalues_between_two_tables(T1,T2)

pvalues = zeros(length(T1.Variables),1);
stats_values = cell(length(T1.Variables),1);
for k = 1: length(T1.Variables)
    d1 = T1(k,:).Variables;
    d2 = T2(k,:).Variables;
    [h,p,ci,stats]=ttest2(d1,d2);
    stats_values{k} = stats;
    pvalues(k) = p;

end

end