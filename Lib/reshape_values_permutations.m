function pValueMatrixRes = reshape_values_permutations(pValueMatrix, combos_all)

N = max(combos_all(:,1));
M = max(combos_all(:,2));
pValueMatrixRes = zeros(N,M);


for i = 1 : size(combos_all,1)
    pValueMatrixRes(combos_all(i,1),combos_all(i,2)) = pValueMatrix(i);
end

end
