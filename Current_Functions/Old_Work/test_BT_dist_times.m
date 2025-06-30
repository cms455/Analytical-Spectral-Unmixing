

num_points = 120;
%A = generate_spectrum_curve(num_points, 3, 200,400,2);
A = build_absorption_matrix(680, 800, [1, 1, 1, 1, 1], num_points);
delta = 14;
[best_indices, frob_norm] = BT_dist_algo_v2(A, 6, 300, delta,2);
disp('Best indices:');
disp(best_indices);
fprintf('Frobenius norm: %.6f\n', frob_norm);
