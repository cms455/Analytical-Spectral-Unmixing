load_A = load('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/combined_spectra_v2.mat');
load_A = load_A.combined_spectra;

load_methyl_blue = readmatrix('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/MethyleneBlue_spectrum.csv');
methyl_blue = load_methyl_blue(3:end,:);
methyl_blue(:,2) = methyl_blue(:,2);


methyl_blue_select = 1e-3*methyl_blue(239:end,2);

load_A_select = load_A(1:121,:);

df = 2;
load_A_select = load_A_select(1:df:end,:);

load_A_select(:,6) = methyl_blue_select(:);

save('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/combined_spectra_w_MB.mat','load_A_select');