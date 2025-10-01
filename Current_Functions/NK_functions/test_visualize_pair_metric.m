min_w = 600;
max_w =1100;
species_bool = [1, 1, 1, 1, 1];
num_points = 150;

k = 5;
wavelengths = linspace(min_w, max_w, num_points);
num_species = sum(species_bool);
%A = build_absorption_matrix(min_w, max_w, species_bool, num_points);
%A = randn(5, 120);  % example matrix
A = generate_spectrum_curve(num_points, 2, min_w, max_w, 5);
pair_map = visualize_pairwise_metric(A);


% -- 3D Surface Plot --

epsilon = 1e-10;
log_metric = log10(pair_map + epsilon);
figure;
[X, Y] = meshgrid(1:num_points, 1:num_points);
surf(X, Y, log_metric, 'EdgeColor', 'none');
xlabel('Column Index i');
ylabel('Column Index j');
zlabel('Frobenius Norm');
title('3D Landscape of Pairwise Frobenius Norms');

view(45, 30);  % Angle the view for better depth