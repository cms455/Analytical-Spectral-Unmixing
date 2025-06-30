function [conditioned_indices, min_inv_indices, submatrix_cond, submatrix_inv, min_cond_val, min_inv_val] = bourgain_tzafriri_k_cols_thresh(A_reg, A,k,num_iter)
    [num_rows, num_cols] = size(A);
    submatrix = [];
    selected_indices = [];
    conditioned_indices = zeros(1,k);
    min_inv_indices = zeros(1,k);
    s=k;
    count = 1;
    inv_val = 2.0;
    cond_val = 2.0;
    while inv_val > 1.1 || count < 10000
        selected_indices = cond_reduce(A, s,num_cols);
        current_cond_val = cond(A_reg(:,selected_indices));
        current_inv_val= norm(pinv(A_reg(:, selected_indices)),'Fro');
        if inv_val > current_inv_val 
            inv_val = current_inv_val;
            min_inv_indices = selected_indices;
        end
        if cond_val > current_cond_val 
            cond_val = current_cond_val;
            conditioned_indices = selected_indices;
        end
        count = count + 1;
    end
    min_cond_val = cond_val;
    min_inv_val = inv_val;
    submatrix_cond = A(:, conditioned_indices);
    submatrix_inv = A(:,min_inv_indices);
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
