data = load("/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/tumor_spectra.mat");

US_data = data.US;
lambda_data = data.PA_lambda(2).image;

figure; imagesc(lambda_data)