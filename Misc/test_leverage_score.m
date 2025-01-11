function selected_columns = leverage_score_column_selection(A, k)
    % A: Input matrix (m x n, e.g., 2 x 150)
    % k: Number of columns to select (e.g., 2)
    
    % Step 1: Compute SVD of A
    [U, S, V] = svd(A); % Economy SVD for efficiency
    
    % Step 2: Compute leverage scores for columns
    r = size(A,2); % Effective rank of A (or use r = min(size(A)))
    V_r = V(:, 1:r); % First r columns of V
    leverage_scores = sum(V_r.^2, 2); % Row-wise squared norms (corresponds to columns of A)
    
    % Step 3: Sort columns by leverage scores
    [~, sorted_indices] = sort(leverage_scores, 'descend');
    
    % Step 4: Select top k columns
    selected_columns = sorted_indices(1:k);
    
    % Optional: Evaluate the norm of the inverse of the selected submatrix
    B = A(:, selected_columns);
    inv_norm = norm(pinv(B), 2); % 2-norm of the inverse
    fprintf('Selected columns: %s\n', mat2str(selected_columns));
    fprintf('Norm of the inverse of the selected submatrix: %f\n', inv_norm);
end

load_A = readmatrix('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/HbOHbR.csv');
A = load_A(3:end, 2:end);
A = A';
shift = 225;
A = A(:, shift:end);
num_rows = size(A,1);
num_cols = size(A,2);
k = 10; % Number of columns to select
selected_columns = leverage_score_column_selection(A, k);
