function [selected_indices, submatrix] = greedy_column_selection(A, Nt)
    % Input:
    % A: Original unstandardized matrix (rows = molecules, cols = wavelengths)
    % Nt: Number of columns to select (size of the subset)
    %
    % Output:
    % selected_indices: Indices of selected columns
    % submatrix: Submatrix formed by the selected columns
    
    [num_rows, num_cols] = size(A);  % Dimensions of the matrix
    selected_indices = [];          % Initialize the list of selected indices
    remaining_indices = 1:num_cols; % Indices of remaining columns
    
    % Iteratively select columns
    for step = 1:Nt
        best_col = -1;              % Initialize the best column index
        best_condition = inf;       % Initialize the best condition number
        
        % Iterate over all remaining columns
        for i = 1:length(remaining_indices)
            current_col = remaining_indices(i);
            trial_indices = [selected_indices, current_col];
            E = A(:, trial_indices);  % Form the submatrix with trial columns
            
            % Compute condition number of the submatrix
            condition_number = cond(E);  % MATLAB's cond() computes cond(E, 2)
            
            % Check if this column improves the condition number
            if condition_number < best_condition
                best_col = current_col;
                best_condition = condition_number;
            end
        end
        
        % Add the best column to the selected set
        selected_indices = [selected_indices, best_col];
        remaining_indices(remaining_indices == best_col) = [];  % Remove it from remaining
    end
    
    % Output the submatrix and selected indices
    submatrix = A(:, selected_indices);
end
