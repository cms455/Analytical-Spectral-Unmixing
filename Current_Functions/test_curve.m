function absorption_curve =  build_curve(num_wavelengths, min_val, max_val, num_peaks, plot_flag)
% Parameters
%num_wavelengths = 100; % Number of wavelengths
%min_wavelength = 400;  % Minimum wavelength (nm)
%max_wavelength = 700;  % Maximum wavelength (nm)
%num_peaks = 3;         % Number of Gaussian peaks

% Generate wavelength range
wavelengths = linspace(min_wavelength, max_wavelength, num_wavelengths);

% Initialize the absorption curve
absorption_curve = zeros(1, num_wavelengths);

% Randomly generate Gaussian peaks
 % For reproducibility
for i = 1:num_peaks
    % Randomly select peak center, height, and width
    peak_center = min_wavelength + (max_wavelength - min_wavelength) * rand();
    peak_height = 0.5 + 1.5 * rand(); % Random height between 0.5 and 2
    peak_width = 10 + 30 * rand();    % Random width between 10 and 40 nm
    
    % Gaussian peak
    peak = peak_height * exp(-((wavelengths - peak_center).^2) / (2 * peak_width^2));
    
    % Add the peak to the absorption curve
    absorption_curve = absorption_curve + peak;
end

% Normalize the absorption curve (optional)
absorption_curve = absorption_curve / max(absorption_curve);
if plot_flag
% Plot the simulated absorption spectrum
figure;
plot(wavelengths, absorption_curve, 'LineWidth', 2);
xlabel('Wavelength (nm)');
ylabel('Absorption');
title('Simulated Absorption Spectrum');
grid on;
end

end