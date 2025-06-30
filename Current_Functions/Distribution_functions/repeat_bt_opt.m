function [best_bt_combo, best_bt_norm] = repeat_bt_opt(A, k, trials1, neigh_range, trials2, N)
% repeat_bt_opt repeats the og_dist_v2 optimization N times
% and returns the best result (lowest bt_norm)
%
% Inputs:
%   A           - input matrix
%   k           - number of columns to select
%   trials1     - number of random search trials
%   neigh_range - neighborhood search radius
%   trials2     - number of neighborhood trials
%   N           - number of optimization repetitions
%
% Outputs:
%   best_bt_combo - best index combination found
%   best_bt_norm  - corresponding minimum norm

    best_bt_norm = Inf;
    best_bt_combo = [];

    for i = 1:N
        [bt_combo, bt_norm] = og_dist_v2(A, k, trials1, neigh_range, trials2);

        if bt_norm < best_bt_norm
            best_bt_norm = bt_norm;
            best_bt_combo = bt_combo;
        end
    end
end


