load_A = readmatrix('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/HbOHbR.csv');
A = load_A(3:end, 2:end);
A = A';
shift = 225;
A = A(:, shift:end);
num_rows = size(A,1);
num_cols = size(A,2);
Nt = 4;
max_cols = 1000;
num_iters = 8;
time_holder=zeros(1,num_iters);
cols  = linspace(2,max_cols,num_iters);
for k= 1:num_iters
k_col = round(cols(k));
A = randn(10,1000);
A_norm = normalize_columns(A);
tic;
[l_submatrix,l_indices] = luke_algorithm(A',2);
time = toc;
time_holder(k) = time;
end


figure;
plot(linspace(2,max_cols,num_iters),time_holder);
xlabel('Number of Iterations');
ylabel('Time')