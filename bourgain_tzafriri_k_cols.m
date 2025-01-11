function [conditioned_indices, submatrix, min_cond_val] = bourgain_tzafriri_k_cols_n_iter(A,k,num_iter)
    [num_rows, num_cols] = size(A);
    submatrix = [];
    ub = floor(log2(num_cols));
    selected_indices = [];
    conditioned_indices = zeros(1,num_cols+1);
    selected_indices_holder = zeros(num_iter, k);
    cond_val_holder = zeros(1, num_iter);
    s=k;
    for m = 1:num_iter
        selected_indices = cond_reduce(A, s,num_cols);
        selected_indices_holder(m,:) = selected_indices
        cond_val_holder(m) = cond(A(:,selected_indices));
    end

    [min_cond_val,min_idx] = min(cond_val_holder);
    conditioned_indices = selected_indices_holder(min_idx,:);
    submatrix = A(:, conditioned_indices);
end




function selected_indices = cond_reduce(A,s,num_cols)
selected_indices = [];
while length(selected_indices) ~= s
random_indices = randperm(num_cols, s);
A_sigma = A(:, random_indices);
G = A_sigma' * A_sigma - eye(s);
alpha = s/4;
F = search_for_F(G,alpha,s);
f = diag(F);
selected_indices = random_indices(f<=(2/s));
end

end

function F = search_for_F(G,alpha,s)

T = floor(alpha^3+1);
J = @(F)[-alpha*F G; G -alpha*F];
J_eig = @(f)eigs([-alpha*diag(f) G; G -alpha*diag(f)],1,'largestreal');
%f = entropic_mirror_descent(J,s,T,alpha);
Aeq = ones(1, s); % Sum(f) = 1
beq = 1;
f0 = (1 / s) * ones(s, 1);
options = optimoptions('fmincon', 'Algorithm', 'interior-point','TolFun',1e-6,'TolX',1e-8);
F_opt = fmincon(J_eig, f0, [], [], Aeq, beq, [], [], [], options);

F = diag(F_opt);
end
