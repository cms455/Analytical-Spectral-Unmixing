function selected_columns = dpp_column_selection(A, k)
    % Step 1: Compute the similarity matrix (kernel)
    L = A' * A; % Use the Gram matrix as the kernel
    
    % Step 2: Eigendecompose the kernel
    [V, D] = eig(L);
    eigenvalues = diag(D);
    
    % Step 3: Sample eigenvectors
    selected_eigenvectors = [];
    for i = 1:length(eigenvalues)
        if rand() < eigenvalues(i) / (1 + eigenvalues(i))
            selected_eigenvectors = [selected_eigenvectors, V(:, i)];
        end
    end
    
    % Step 4: Greedy column selection
    selected_columns = [];
    for i = 1:k
        scores = sum((A' .* selected_eigenvectors).^2, 2);
        [~, best_col] = max(scores);
        selected_columns = [selected_columns, best_col];
        selected_eigenvectors = selected_eigenvectors - ...
            (selected_eigenvectors .* A(:, best_col)') .* A(:, best_col)';
    end
    
    % Ensure uniqueness
    selected_columns = unique(selected_columns);
    selected_columns = selected_columns(1:k); % Limit to k columns
end

% Example Usage
A = rand(2, 150); % Example matrix
k = 2; % Number of columns to select
selected_columns = dpp_column_selection(A, k);
disp('Selected columns:');
disp(selected_columns);
