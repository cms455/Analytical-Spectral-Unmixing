function [bestIndices, minNorm] = brute_force(A, k)
    [m, n] = size(A);

    if k > n
        error('k must be less than or equal to the number of columns in A');
    end

    combos = nchoosek(1:n, k);
    numCombos = size(combos, 1);

    minNorm = inf;
    bestIndices = [];

    for i = 1:numCombos
        cols = combos(i, :);
        submatrix = A(:, cols);
        pinvSub = pinv(submatrix);
        normVal = norm(pinvSub, 'fro');

        if normVal < minNorm
            minNorm = normVal;
            bestIndices = cols;
        end
    end
end
