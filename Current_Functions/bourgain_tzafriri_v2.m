function [conditioned_indices, submatrix] = bourgain_tzafriri_v2(A)
    [num_rows, num_cols] = size(A);  
    submatrix = [];
    ub = floor(log2(num_cols));
    selected_indices = [];
    conditioned_indices = zeros(1,num_cols+1);
    for si = 1:ub
        s = 2^si;
        for k = 1:8*log2(s)
            selected_indices = cond_reduce(A, s,num_cols);
               
            if cond(A(:,selected_indices)) <= sqrt(3)-0.6
                conditioned_indices = selected_indices;
                break
            end
            
        end
        if length(conditioned_indices) < s
                submatrix = A(:, conditioned_indices);
                return
        end
    end
   
end




function selected_indices = cond_reduce(A,s,num_cols)
random_indices = randperm(num_cols, s);
A_sigma = A(:, random_indices);
G = A_sigma' * A_sigma - eye(s);
alpha = s/4;
F = search_for_F(G,alpha,s);
f = diag(F);
selected_indices = random_indices(f<=(2/s));
end

function F = search_for_F(G,alpha,s)

T = floor(alpha^3+1);
J = @(F)[-alpha*F G; G -alpha*F];
J_eig = @(f)eigs([-alpha*diag(f) G; G -alpha*diag(f)],1,'largestreal');
f = entropic_mirror_descent(J,s,T,alpha);
F = diag(f);
%{
Aeq = ones(1, s); % Sum(f) = 1
beq = 1;
f0 = (1 / s) * ones(s, 1);
options = optimoptions('fmincon', 'Algorithm', 'interior-point','TolFun',1e-6,'TolX',1e-8);
F_opt = fmincon(J_eig, f0, [], [], Aeq, beq, [], [], [], options);
F = diag(F_opt);
%}
end
