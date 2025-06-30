% Main script to run pandemicInstance repeatedly until epidemic is large

R0 = 0.9999;
I0 = 1;
threshold = 1e5;

while true
    [I_vec, S] = pandemicInstance(R0, I0);
    if S >= threshold
        break;
    end
end

% Plot the epidemic trajectory
figure;
plot(0:length(I_vec)-1, I_vec, '-o', 'LineWidth', 1.5);
xlabel('Generation n');
ylabel('Number Infected I_n');
title(sprintf('Pandemic Instance (R_0 = %.4f, Total Infected = %.0f)', R0, S));
grid on;
