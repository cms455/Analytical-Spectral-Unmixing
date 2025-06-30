function A = build_absorption_matrix_v2(min_w, max_w, species_bool, num_points)

load_A = load('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/combined_spectra_v2.mat');
load_A = load_A.combined_spectra;
wavelengths = load_A(:,1);
full_A = load_A(:,:)';

%{
HbO_data = full_A(1,:);
Hb_data = full_A(2,:);
ICG_data = full_A(3,:);
lipid_data = full_A(4,:);

%}

HbO_data = full_A([1,2],:);
Hb_data = full_A([1,3],:);
ICG_data = full_A([1,4],:);
lipid_data = full_A([1,5],:);

load_methyl_blue = readmatrix('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/MethyleneBlue_spectrum.csv');
methyl_blue = load_methyl_blue(3:end,:);
methyl_blue(:,2) = 0.1*methyl_blue(:,2);

species_data = {Hb_data, HbO_data, lipid_data, ICG_data, methyl_blue};

num_species = sum(species_bool);

A = zeros(num_species,num_points);
% species_bool is 1 if you want to include the species, and 0 if you want
% to skip.
% Order, [Hb, HbO, ICG, Lipid]


count = 1;
for n = 1:length(species_bool)
    if species_bool(n)
        curve = species_data{n};
        selected_curve_interp = pick_bandwidth(curve, min_w,max_w, num_points);
        A(count, :) = selected_curve_interp;
        count = count+ 1;
    end

   
end


end





function selected_curve_interp = pick_bandwidth(curve, min_w, max_w, num_points)
start_idx = find(curve(:,1) == min_w );
end_idx = find(curve(:,1) == max_w);
wavelengths = start_idx:end_idx;
dense_wavelengths = linspace(start_idx,end_idx, num_points);
selected_curve = curve(start_idx:end_idx,2);
selected_curve_interp = interp1(wavelengths, selected_curve, dense_wavelengths,'Spline');
%{
figure; 
hold on;
plot(wavelengths,selected_curve,'ro');
plot(dense_wavelengths,selected_curve_interp,'bo');
%}
end

