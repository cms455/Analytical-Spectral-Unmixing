function S_list = pandemicEnsemble(R0, I0, N)
% Runs N simulations and returns a vector of epidemic sizes
% INPUTS:
%   R0  - reproduction number
%   I0  - initial infections
%   N   - number of simulations
% OUTPUT:
%   S_list - array of total sizes S for each epidemic

    S_list = zeros(1, N);
    for i = 1:N
        S_list(i) = runOne(R0, I0);
    end
end

function S = runOne(R0, I0)
% Helper to compute a single epidemic size
    max_generations = 1e4;
    I = I0;
    S = I0;

    for n = 2:max_generations
        mean_next = R0 * I;
        I = poissrnd(mean_next);
        S = S + I;
        if I == 0
            break;
        end
    end
end
