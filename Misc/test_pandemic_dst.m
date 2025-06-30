function epidemicHistograms()
    % Parameters
    R0 = 0.99;
    I0 = 1;
    N = 1e4;
    num_log_bins = 50;

    % Generate ensemble
    S_list = pandemicEnsemble(R0, I0, N);
    S_list = S_list(S_list > 0);  % remove zeros

    %% --- Part (e): Linear-binned histogram ---
    figure;
    subplot(1,2,1);
    histogram(S_list, 100, 'Normalization', 'pdf');
    xlabel('Epidemic size S');
    ylabel('P(S)');
    title(sprintf('Linear-Binned Histogram (R_0 = %.2f)', R0));
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    grid on;

    %% --- Part (f): Logarithmic-binned histogram ---
    subplot(1,2,2);

    % Log-binning
    min_S = min(S_list);
    max_S = max(S_list);
    edges = unique(round(logspace(log10(min_S), log10(max_S), num_log_bins)));
    if numel(edges) < 2
        error('Too few unique bins. Adjust binning or epidemic size.');
    end

    bin_counts = histcounts(S_list, edges);
    bin_centers = sqrt(edges(1:end-1) .* edges(2:end));
    bin_widths = diff(edges);
    P_S = bin_counts ./ (bin_widths * numel(S_list));

    % Plot log-binned data
    loglog(bin_centers, P_S, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
    hold on;

    % Power-law overlay
    tau = 3/2;
    S_fit = logspace(log10(min(bin_centers)), log10(max(bin_centers)), 200);
    C = P_S(1) * bin_centers(1)^tau;
    loglog(S_fit, C * S_fit.^(-tau), 'r--', 'LineWidth', 2);

    xlabel('Epidemic size S');
    ylabel('P(S)');
    title('Log-Binned Histogram + Power Law');
    legend('Simulated data', '\tau = 3/2', 'Location', 'southwest');
    grid on;
        %% --- Part (g): Scale by S^{3/2} to test power-law collapse ---
    figure;
    scaled_P = P_S .* (bin_centers .^ (3/2));  % scale P(S) by S^{3/2}

    loglog(bin_centers, scaled_P, 'ko-', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'LineWidth', 1.5);
    xlabel('Epidemic size S');
    ylabel('Scaled distribution: S^{3/2} P(S)');
    title('Scaling Collapse: S^{3/2} P(S) vs. S');
    grid on;

end

function S_list = pandemicEnsemble(R0, I0, N)
    % Runs N simulations and returns a vector of total epidemic sizes
    S_list = zeros(1, N);
    for i = 1:N
        S_list(i) = runOne(R0, I0);
    end
end

function S = runOne(R0, I0)
    % Simulates a single epidemic trajectory and returns total infected S
    max_generations = 1e4;
    I = I0;
    S = I0;

    for n = 2:max_generations
        I = poissrnd(R0 * I);
        S = S + I;
        if I == 0
            break;
        end
    end
end
