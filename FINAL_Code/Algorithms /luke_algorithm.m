function [selected_matrix,row_indices] = luke_algorithm(A, k)
 
    [n, ~] = size(A);

    if k > n
        error('k must be less than or equal to the number of rows in the matrix.');
    end
    

    row_indices = 1:n; 
    
   
    while length(row_indices) > k
        best_sigma_min = -inf;  
        worst_row_to_remove = -1;  %
        
        for i = 1:length(row_indices)
         
            temp_indices = row_indices;
            temp_indices(i) = [];
            temp_matrix = A(temp_indices, :);
     
            sigma_min = min(svd(temp_matrix));
            
           
            if sigma_min > best_sigma_min
                best_sigma_min = sigma_min;
                worst_row_to_remove = i;
            end
        end
       
        row_indices(worst_row_to_remove) = [];
    end
   
    selected_matrix = A(row_indices, :);
end
