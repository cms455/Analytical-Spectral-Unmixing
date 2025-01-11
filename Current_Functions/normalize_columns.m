function standardized_matrix = normalize_columns(matrix)
    % Input:
    % matrix: A numerical matrix (n x m)
    %
    % Output:
    % standardized_matrix: Matrix with each column normalized to unit l2 norm

    % Compute the l2-norm of each column
    column_norms = sqrt(sum(matrix.^2, 1));  % 1 x m vector of column norms

    % Avoid division by zero (if a column is all zeros, its norm will be zero)
    column_norms(column_norms == 0) = 1;

    % Normalize each column
    standardized_matrix = matrix ./ column_norms;
end