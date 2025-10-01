function pairwise_metric = visualize_pairwise_metric(A)
% visualize_pairwise_metric: Plots metric (e.g., Frobenius norm of pseudoinverse)
% for all column pairs (i,j) of matrix A.
%
% Input:
%   A - (m x n) input matrix
%
% Output:
%   pairwise_metric - (n x n) matrix of metric values for all pairs

[m, n] = size(A);
pairwise_metric = nan(n, n);  % Initialize with NaNs to avoid plotting diagonal

for i = 1:n
    for j = i+1:n
        idx = [i, j];
        submatrix = A(:, idx);
        value = norm(pinv(submatrix), 'fro');
        pairwise_metric(i, j) = value;
        pairwise_metric(j, i) = value;  % symmetry
    end
end

% Plot heatmap with log-scaled colorbar
figure;

% Take log10 of the matrix, adding a small epsilon to avoid log(0)
epsilon = 1e-10;
log_metric = log10(pairwise_metric + epsilon);
imagesc(log_metric);
axis square;

% Adjust colorbar ticks to reflect the original scale
colormap(jet);  % or colormap(parula) for a perceptually better scale
c = colorbar();
caxis([min(log_metric(:)), max(log_metric(:))]);


% Choose tick values in log space, then label with real values
log_ticks = floor(min(log_metric(:))):ceil(max(log_metric(:)));
c.Ticks = log_ticks;
c.TickLabels = arrayfun(@(x) sprintf('%.1e', 10.^x), log_ticks, 'UniformOutput', false);

xlabel('Column Index');
ylabel('Column Index');
title('Pairwise Frobenius Norm of Pseudoinverse (Log Scale)');
