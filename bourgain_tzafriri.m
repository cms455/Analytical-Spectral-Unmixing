function [selected_indices, submatrix] = bourgain_tzafriri(A, Nt, max_iter)
    % Input:
    % A: Standardized input matrix (rows = molecules, cols = wavelengths)
    % Nt: Number of columns to select
    % max_iter: Maximum number of iterations for refinement
    
    [num_rows, num_cols] = size(A);  % Dimensions of the matrix
    if Nt > num_cols
        error('Nt must be less than or equal to the number of columns in A');
    end
    

    selected_indices = [];
    best_condition = inf; 
    max_s = min(2 * Nt, num_cols);  

    for s = 4:4:max_s  
        for iter = 1:max_iter
            random_indices = randperm(num_cols, s);
            A_sigma = A(:, random_indices);  
            
            G = A_sigma' * A_sigma - eye(s);  % Hollow Gram matrix
            
            [U, D, ~] = svd(G); 
            diag_d = diag(D);   
            
            refined_indices = random_indices(diag_d.^2 <= (2));
            
            if length(refined_indices) >= Nt

                A_tau = A(:, refined_indices(1:Nt));
                current_condition = cond(A_tau);
               
                if current_condition < best_condition
                    best_condition = current_condition;
                    selected_indices = refined_indices(1:Nt);
                end
            end
        end
    end
    submatrix = A(:, selected_indices);
end
