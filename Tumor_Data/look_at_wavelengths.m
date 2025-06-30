min_w = 500;
max_w = 600;
species_bool = [1, 1, 0, 0, 0];
num_points = 100;
wavelengths = linspace(min_w, max_w, num_points);
num_species = sum(species_bool);
A = build_absorption_matrix(min_w, max_w, species_bool, num_points);

num_points_peek = 40;

hb_val_holder = zeros(num_points_peek,1);

for i = 1:num_points_peek
    hb_val_holder(i) = cond(A(:,[60,i]));
end

figure();
plot(wavelengths(1:num_points_peek),hb_val_holder);