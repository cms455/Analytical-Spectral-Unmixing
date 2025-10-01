function [best_genotype, best_fitness] = nk_column_selector_anneal(A, k, num_iters, T_init, T_min, alpha, kB)
% nk_column_selector_anneal: Simulated annealing for column selection
%
% Inputs:
%   A          - Input matrix (m x n)
%   k          - Number of columns to select
%   num_iters  - Total number of iterations
%   T_init     - Initial temperature
%   T_min      - Minimum temperature (stopping condition)
%   alpha      - Cooling rate (e.g. 0.95)
%   kB         - Boltzmann constant (e.g. 1)
%
% Outputs:
%   best_genotype - Binary vector indicating selected columns
%   best_fitness  - Negative Frobenius norm

n = size(A, 2);

% Initial genotype with k random 1s
genotype = zeros(1, n);
rand_idx = randperm(n, k);
genotype(rand_idx) = 1;

% Initial fitness
fitness = -norm(pinv(A(:, genotype == 1)), 'fro');

% Store best
best_genotype = genotype;
best_fitness = fitness;

% Initialize temperature
T = T_init;

for iter = 1:num_iters
    % Swap one 1 with one 0
    ones_idx = find(genotype == 1);
    zeros_idx = find(genotype == 0);
    i = ones_idx(randi(length(ones_idx)));
    j = zeros_idx(randi(length(zeros_idx)));
    
    new_genotype = genotype;
    new_genotype(i) = 0;
    new_genotype(j) = 1;
    
    % Compute new fitness
    selected = find(new_genotype);
    new_fitness = -norm(pinv(A(:, selected)), 'fro');
    
    delta_E = new_fitness - fitness;

    if delta_E > 0
        % Accept improvement
        genotype = new_genotype;
        fitness = new_fitness;
    else
        % Accept with probability exp(delta_E / (kB * T))
        if rand() < exp(delta_E / (kB * T))
            genotype = new_genotype;
            fitness = new_fitness;
        end
    end
    
    % Update best if needed
    if fitness > best_fitness
        best_genotype = genotype;
        best_fitness = fitness;
    end
    
    % Cool down
    T = max(T_min, alpha * T);
end
end
