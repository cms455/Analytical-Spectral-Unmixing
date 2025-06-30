load_A = readmatrix('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/HbOHbR_Spectrum.csv');
A = load_A(3:end, 2:end);
A = A';
shift = 225;
A = A(:, shift:end);
num_rows = size(A,1);
num_cols = size(A,2);
Nt = 4;
max_iter = 100;

k_items = 2;
x_idx = 1:length(A);
combinations = nchoosek(x_idx,k_items);

tic;
holder = zeros(1,length(combinations));
for i = 1:length(combinations)
    idx = combinations(i,:);
    holder(i) = norm(pinv(A(:,idx)),'Fro');
end
toc;
[min_val,min_idx] = min(holder,[],'all');
disp(combinations(min_idx,:));
disp('Inverse Val:')
disp(holder(min_idx));