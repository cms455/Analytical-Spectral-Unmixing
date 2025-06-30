
% Generate wavelength range
min_wavelength = 400; % Example value, update to match your input
max_wavelength = 700; % Example value, update to match your input
num_species = 3;
num_wavelengths = 300;
wavelengths = linspace(min_wavelength, max_wavelength, num_wavelengths);

A = generate_spectrum_curve(num_wavelengths, num_species, 400, 700, 3);

% Plot each row of A
figure;
hold on;
for i = 1:num_species
    plot(wavelengths, A(i, :), 'LineWidth', 2);
end
hold off;

% Customize the plot
xlabel('Wavelength (nm)');
ylabel('Absorption');
title('Absorption Spectra');
legend(arrayfun(@(x) sprintf('Species %d', x), 1:num_species, 'UniformOutput', false));
grid on;
