function [muvalues,md1values,md2values] = compute_meanstd_between_two_tables(T1,T2)

muvalues = zeros(length(T1.Variables),1);
md1values = zeros(length(T1.Variables),1);
md2values = zeros(length(T1.Variables),1);
for k = 1: length(T1.Variables)
    d1 = T1(k,:).Variables;
    d2 = T2(k,:).Variables;
    md1= mean(d1);
    md2 = mean(d2);
    mu = mean(d1)-mean(d2);
    muvalues(k) = mu;
    md1values(k) = md1;
    md2values(k) = md2;
end

end