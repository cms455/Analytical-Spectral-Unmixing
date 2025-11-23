function [best_genotype, best_fitness] = nk_column_selector(A, k, num_iters)
% nk_column_selector: NK-model style local search for minimizing Frobenius norm
%
% Inputs:
%   A          - Input matrix (m x n)
%   k          - Number of columns to select (number of 1s in genotype)
%   num_iters  - Number of local search iterations
%
% Outputs:
%   best_genotype - Binary vector (1 x n), k ones indicate selected columns
%   best_fitness  - Negative Frobenius norm (higher is better fitness)

n = size(A, 2);  % total number of columns

% Initialize with k random 1s
genotype = zeros(1, n);
rand_idx = randperm(n, k);
genotype(rand_idx) = 1;

% Evaluate fitness
fitness = -norm(pinv(A(:, genotype == 1)), 'fro');  % negative norm = fitness

% save best
best_genotype = genotype;
best_fitness = fitness;

% local walk (mutation: swap a 1 and a 0)
for iter = 1:num_iters
    % propose mutation: swap one 1 with one 0
    ones_idx = find(genotype == 1);
    zeros_idx = find(genotype == 0);
    
    i = ones_idx(randi(length(ones_idx)));
    j = zeros_idx(randi(length(zeros_idx)));
    
    new_genotype = genotype;
    new_genotype(i) = 0;
    new_genotype(j) = 1;
    
    % evaluate fitness
    selected = find(new_genotype);
    new_fitness = -norm(pinv(A(:, selected)), 'fro');
    
    % change if better
    if new_fitness > best_fitness
        best_genotype = new_genotype;
        best_fitness = new_fitness;
        genotype = new_genotype;  % move to new point
    end
end
%{
selected_indices = find(best_genotype);
fprintf('Selected column indices: [%s]\n', num2str(selected_indices));
fprintf('Best Fitness (−Frobenius norm): %.8f\n', best_fitness);
% Visualization
figure;
imagesc(A(:, selected_indices));
xlabel('Selected Columns'); ylabel('Row Index');
title(sprintf('Submatrix A(:, selected) — Fitness = %.8f', best_fitness));
colorbar;
%}
end
