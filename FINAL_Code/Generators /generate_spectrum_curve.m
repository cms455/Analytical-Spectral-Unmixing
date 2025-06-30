function A = generate_spectrum_curve(num_wavelengths, num_species, min_wavelength, max_wavelength, num_peaks)
    plot_flag = false;
    A = zeros(num_species, num_wavelengths);
    for i = 1:num_species
        absorption_curve = build_curve(num_wavelengths, min_wavelength,max_wavelength, num_peaks, plot_flag);
        A(i,:) = absorption_curve;
    end


end