function f = entropic_mirror_descent_v2(J, s, T, alpha)
    f_1 = ones(s, 1) / s;
    f_t = repmat(f_1, 1, T + 1);
    eig_val_holder = Inf(1, T);

    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('local'); % ensure pool exists
    end

    for t = 1:T
        try
            fut = parfeval(@eigs, 2, J(diag(f_t(:, t))), 1, 'largestreal');
            startTime = tic;
            while toc(startTime) < 1
                if strcmp(fut.State, 'finished')
                    [V, D] = fetchOutputs(fut);
                    eig_val_holder(t) = D;
                    break;
                end
                pause(0.01); % allow time for future to complete
            end

            if strcmp(fut.State, 'finished') == false
                cancel(fut);
                warning('Iteration %d: eigs timed out. Skipping.', t);
                f_t(:, t + 1) = f_t(:, t);
                continue;
            end

            u = V(1:length(V)/2);
            v = V(length(V)/2+1:end);
            theta = -alpha * (u.^2 + v.^2);

            if norm(theta, 'inf') == 0
                f_t(:, t + 1) = f_t(:, t);
                continue;
            end

            beta = sqrt(2 * log(s) / (T * norm(theta, 'inf')^2));
            h = f_t(:, t) .* exp(-beta .* theta);
            f_t(:, t + 1) = h / sum(h);

        catch ME
            warning('Iteration %d failed: %s', t, ME.message);
            f_t(:, t + 1) = f_t(:, t);  % carry forward on error
        end
    end

    [~, min_idx] = min(eig_val_holder);
    f = f_t(:, min_idx);
end
