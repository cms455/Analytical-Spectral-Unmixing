function [selected_matrix,row_indices] = luke_algorithm(matrix, num_rows_to_select)
    % Input:
    % matrix: The full matrix with dimensions (n x m), where n is the number of rows
    % num_rows_to_select: The number of rows to retain
    %
    % Output:
    % selected_matrix: The resulting submatrix with the selected rows
    
    % Get the total number of rows in the matrix
    [n, ~] = size(matrix);
    
    % Validate inputs
    if num_rows_to_select > n
        error('num_rows_to_select must be less than or equal to the number of rows in the matrix.');
    end
    
    % Initialize the working set of rows
    row_indices = 1:n;  % Start with all rows
    
    % Iteratively remove rows until the desired number of rows is reached
    while length(row_indices) > num_rows_to_select
        best_sigma_min = -inf;  % Initialize the best sigma_min found
        worst_row_to_remove = -1;  % Initialize the worst row index to remove
        
        % Try removing each row and compute the smallest singular value
        for i = 1:length(row_indices)
            % Create a temporary matrix by excluding the i-th row
            temp_indices = row_indices;
            temp_indices(i) = [];
            temp_matrix = matrix(temp_indices, :);
            
            % Compute the smallest singular value of the temporary matrix
            sigma_min = min(svd(temp_matrix));
            
            % Check if this is the best row to remove
            if sigma_min > best_sigma_min
                best_sigma_min = sigma_min;
                worst_row_to_remove = i;
            end
        end
        
        % Remove the worst row from the working set
        row_indices(worst_row_to_remove) = [];
    end
    
    % Construct the resulting submatrix
    selected_matrix = matrix(row_indices, :);
end
