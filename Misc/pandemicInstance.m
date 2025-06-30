function [I_vec, S] = pandemicInstance(R0, I0)
% Simulates a single pandemic instance.
% INPUTS:
%   R0  - basic reproduction number
%   I0  - initial number of infected individuals
% OUTPUTS:
%   I_vec - vector of infected counts per generation
%   S     - total number of infections over all generations

    max_generations = 1e4;
    I_vec = zeros(1, max_generations);
    I_vec(1) = I0;

    for n = 2:max_generations
        mean_next = R0 * I_vec(n-1);
        I_vec(n) = poissrnd(mean_next);
        if I_vec(n) == 0
            I_vec = I_vec(1:n);
            break;
        end
    end

    S = sum(I_vec);
end
