min_w = 680;
max_w =1000;
species_bool = [1, 1, 1, 0, 0];
num_points = 150;
wavelengths = linspace(min_w, max_w, num_points);
num_species = sum(species_bool);
A = build_absorption_matrix(min_w, max_w, species_bool, num_points);

delta = 2;
for i = 1:10
[conditioned_indices, submatrix] = bourgain_tzafriri_v2(A,delta);
k_A = cond(A);
nK_A = cond(A(:,conditioned_indices));
nE_A = norm(pinv(A(:,conditioned_indices)),'Fro');
fprintf("\n OG K: %d \n",k_A);
fprintf("New K: %d ",nK_A);
fprintf("\nNew E: %d",nE_A);
disp("Indices ")
disp(conditioned_indices);



end
