% Generate wavelength range
min_wavelength = 400; % Example value
max_wavelength = 700; % Example value
num_species = 2;
num_wavelengths = 50;
wavelengths = linspace(min_wavelength, max_wavelength, num_wavelengths);

% Generate the spectrum matrix A
A = generate_spectrum_curve(num_wavelengths, num_species, min_wavelength, max_wavelength, 3);
A_norm = normalize_columns(A);
k_items = 4;
x_idx = 1:length(A);
num_iters = 5000;
brute_search_val_holder = zeros(k_items -1);

for k = 2:k_items
    combinations = nchoosek(x_idx, k);

    disp('BF TIME:')
    tic;
    holder = zeros(1, length(combinations));
    for i = 1:length(combinations)
        idx = combinations(i, :);
        holder(i) = norm(pinv(A(:, idx)), 'Fro');
    end
    toc;


    [min_val, min_idx] = min(holder, [], 'all');
    selected_indices = combinations(min_idx, :);
    disp('Selected Indices BF:');
    disp(selected_indices);
    disp('Inverse Val BF:');
    disp(holder(min_idx));
    brute_search_val_holder(k-1) = holder(min_idx);
    
    disp('BT ALGO')
    tic;
    [conditioned_indices, min_inv_indices, submatrix_cond, submatrix_inv, min_cond_val, min_inv_val] = bourgain_tzafriri_all_fix_selections(A, A_norm, k, num_iters);
    toc;

    disp('Selected Indices BT:')
    disp(min_inv_indices)
    disp('Inverse Val BT:')
    disp(min_inv_val)


    disp('Luke ALGO')
    tic; 
    [l_submatrix,l_indices] = luke_algorithm(A',k);
    toc;


    disp('Selected Indices LUKE:')
    disp(l_indices)
    disp('Inverse Val LUKE:')
    min_luke_val = norm(pinv(l_submatrix),'Fro');
    disp(min_luke_val)




    % Plot all rows of A
    figure;
    hold on;
    for i = 1:num_species
        plot(wavelengths, A(i, :), 'LineWidth', 2); % Plot all rows
    end
    
    % Add red lines at selected indices
    for i = 1:length(selected_indices)
        x_value = wavelengths(selected_indices(i));
        xline(x_value, 'r--', 'LineWidth', 1.5); % Use yline for vertical lines
        x_value = wavelengths(min_inv_indices(i));
        xline(x_value, 'b--', 'LineWidth', 1.5); % Use yline for vertical lines
    end
    hold off;

    % Customize the plot
    xlabel('Wavelength (nm)');
    ylabel('Absorption');
    title('All Absorption Spectra with Selected Indices Highlighted');
    legend(arrayfun(@(x) sprintf('Species %d', x), 1:num_species, 'UniformOutput', false), 'Location', 'Best');
    grid on;
end
