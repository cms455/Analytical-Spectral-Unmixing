load_A = readmatrix('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/HbOHbR.csv');
A = load_A(3:end, 2:end);
A = A';
shift = 225;
A = A(:, shift:end);
num_rows = size(A,1);
num_cols = size(A,2);
Nt = 4;


tic;
[l_submatrix,l_indices] = luke_algorithm(A',4);
toc;

%disp(l_submatrix)
disp(l_indices);

test = A(:,l_indices);
fprintf('Condition Number:%d \n', cond(test));
fprintf('Inverse norm: %d \n', norm(pinv(A(:,l_indices)),'Fro'));