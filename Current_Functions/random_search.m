function [selected_indices, min_inv_val] = random_search(A,k, num_iters)

holder_indices = zeros(num_iters,k);
holder_min_inv = zeros(1,num_iters);
num_cols = size(A,2);

for i = 1:num_iters
 indices = squeeze(randi(num_cols,1,k));
 while length(unique(indices)) ~= k
     indices = squeeze(randi(num_cols,1,k)); 
 end
 submatrix = A(:,indices);
 norm_inv = norm(pinv(submatrix),'Fro');
 holder_min_inv(i) = norm_inv;
 holder_indices(i,:) = indices;
end

[min_val,min_idx] = min(holder_min_inv);

selected_indices = holder_indices(min_idx,:);
min_inv_val = min_val;
end

