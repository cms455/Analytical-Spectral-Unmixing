function [conditioned_indices, min_inv_indices, submatrix_cond, submatrix_inv, min_cond_val, min_inv_val] = bourgain_tzafriri_all_fix_selections(A_reg, A,k,num_iter)
    [num_rows, num_cols] = size(A);
    submatrix = [];
    ub = floor(log2(num_cols));
    selected_indices = [];
    conditioned_indices = zeros(1,num_cols+1);
    selected_indices_holder = zeros(num_iter, k);
    cond_val_holder = zeros(1, num_iter);
    inv_val_holder = zeros(1, num_iter);
    s=k;
    for m = 1:num_iter
        selected_indices = cond_reduce(A, s,num_cols);
        selected_indices_holder(m,:) = selected_indices;
        cond_val_holder(m) = cond(A_reg(:,selected_indices));
        inv_val_holder(m) = norm(pinv(A_reg(:, selected_indices)),'Fro');
    end

    [min_cond_val,min_cond_idx] = min(cond_val_holder);
    [min_inv_val,min_inv_idx] = min(inv_val_holder);
    conditioned_indices = selected_indices_holder(min_cond_idx,:);
    min_inv_indices = selected_indices_holder(min_inv_idx,:);
    submatrix_cond = A(:, conditioned_indices);
    submatrix_inv = A(:,min_inv_indices);
end




function selected_indices = cond_reduce(A,s,num_cols)
selected_indices = [];
empty_count = 0;
s_original = s;
while length(selected_indices) ~= s_original
if(empty_count > 2)
    s=s+1;
    empty_count = 0;
    %disp(s)
end
random_indices = randperm(num_cols, s);
A_sigma = A(:, random_indices);
G = A_sigma' * A_sigma - eye(s);
alpha = s/4;
F = search_for_F(G,alpha,s);
f = diag(F);
selected_indices = random_indices(f<=(6/s));
if length(selected_indices) < s_original
    empty_count = empty_count + 1;
end
end

end

function F = search_for_F(G,alpha,s)

T = floor(alpha^2+1);
J = @(F)[-alpha*F G; G -alpha*F];
J_eig = @(f)eigs([-alpha*diag(f) G; G -alpha*diag(f)],1,'largestreal');
f = entropic_mirror_descent(J,s,T,alpha);
F = diag(f);

end
