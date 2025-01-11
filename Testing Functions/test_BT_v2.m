load_A = readmatrix('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/HbOHbR.csv');
A = load_A(3:end, 2:end);
A = A';
shift = 225;
A = A(:, shift:end);
num_rows = size(A,1);
num_cols = size(A,2);
A = randn(10,50);
A_norm = normalize_columns(A);
tic;
[selected_indices,submatrix] = bourgain_tzafriri_v2(A_norm);
toc;

disp(selected_indices);
A_norm_selected = A_norm(:,selected_indices);
A_selected = A(:,selected_indices);
fprintf('Normalized Condition Number %d \n',cond(A_norm_selected));
fprintf('Condition Number %d \n',cond(A_selected));
fprintf('Inverse norm: %d \n', norm(pinv(A_selected)));
