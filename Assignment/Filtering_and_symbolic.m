% Does the corrected filtering of FFT and also plots the symbolic function. Needs to be
% run after the 'mean_plots.m' file.
% Developed by Aswath Ashok forthe Data Analysis Assignment


% Calculate the FFT of the functions
fft_u = fft(u_subset);

% Find peak frequency for each subplot
[~, peak_index_original] = max(abs(fft_u));
peak_frequency_original = frequencies(peak_index_original);

% Moving Average filter
windowSize = 100;  % Adjust as needed
smoothed_fft_u_movmean = smoothdata(abs(fft_u), 'movmean', windowSize);
[~, peak_index_movmean] = max(smoothed_fft_u_movmean);
peak_frequency_movmean = frequencies(peak_index_movmean);

% Chebyshev filter
rp = 3;  % Passband ripple in dB
rs = 60; % Stopband attenuation in dB
fc = 0.1; % Cutoff frequency in normalized units (0 to 1)
[b, a] = cheby1(6, rp, fc);
smoothed_fft_u_cheby = filter(b, a, abs(fft_u));
[~, peak_index_cheby] = max(smoothed_fft_u_cheby);
peak_frequency_cheby = frequencies(peak_index_cheby);

% Apply Hann window directly to the time-domain signal
u_subset_hann = u_subset .* hann(length(u_subset));
smoothed_fft_u_hann = fft(u_subset_hann);
[~, peak_index_hann] = max(abs(smoothed_fft_u_hann));
peak_frequency_hann = frequencies(peak_index_hann);

% Plot the results in subplots
figure;

% Original FFT
subplot(4, 1, 1);
plot(frequencies, abs(fft_u(1:length(frequencies))), 'b-', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Original FFT');
xlim([0 5]);

% Smoothed FFT using Moving Average filter
subplot(4, 1, 2);
plot(frequencies, smoothed_fft_u_movmean(1:length(frequencies)), 'r-', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Smoothed FFT (Moving Average)');
xlim([0 5]);

% Smoothed FFT using Chebyshev filter
subplot(4, 1, 3);
plot(frequencies, smoothed_fft_u_cheby(1:length(frequencies)), 'g-', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Smoothed FFT (Chebyshev)');
xlim([0 5]);

% Include the Hann window result
subplot(4, 1, 4);
plot(frequencies, abs(smoothed_fft_u_hann(1:length(frequencies))), 'm-', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Smoothed FFT (Hann Window)');
xlim([0 5]);

sgtitle('Comparison of Different Smoothing Techniques for FFT of u');


% Plot all results in subplots with a line of slope -5/3
figure;

% Original FFT
subplot(2, 2, 1);
loglog(frequencies, abs(fft_u(1:length(frequencies))), 'b-', 'LineWidth', 1.5);
hold on;
loglog(frequencies, frequencies.^(-5/3), 'k--', 'LineWidth', 1.5);
scatter(peak_frequency_original, abs(fft_u(peak_index_original)), 'r', 'filled');
text(peak_frequency_original, abs(fft_u(peak_index_original)), ...
    sprintf('(%.2f Hz, %.2f Magnitude)', peak_frequency_original, abs(fft_u(peak_index_original))), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'r');
xlabel('log of Frequency (Hz)');
ylabel('log of Magnitude');
title('Original FFT');
grid on;
hold off;

% Smoothed FFT using Moving Average filter
subplot(2, 2, 2);
loglog(frequencies, smoothed_fft_u_movmean(1:length(frequencies)), 'r-', 'LineWidth', 1.5);
hold on;
loglog(frequencies, frequencies.^(-5/3), 'k--', 'LineWidth', 1.5);
scatter(peak_frequency_movmean, smoothed_fft_u_movmean(peak_index_movmean), 'r', 'filled');
text(peak_frequency_movmean, smoothed_fft_u_movmean(peak_index_movmean), ...
    sprintf('(%.2f Hz, %.2f Magnitude)', peak_frequency_movmean, smoothed_fft_u_movmean(peak_index_movmean)), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'r');
xlabel('log of Frequency (Hz)');
ylabel('log of Magnitude');
title('Smoothed FFT (Moving Average)');
grid on;
hold off;

% Smoothed FFT using Chebyshev filter
subplot(2, 2, 3);
loglog(frequencies, smoothed_fft_u_cheby(1:length(frequencies)), 'g-', 'LineWidth', 1.5);
hold on;
loglog(frequencies, frequencies.^(-5/3), 'k--', 'LineWidth', 1.5);
scatter(peak_frequency_cheby, smoothed_fft_u_cheby(peak_index_cheby), 'r', 'filled');
text(peak_frequency_cheby, smoothed_fft_u_cheby(peak_index_cheby), ...
    sprintf('(%.2f Hz, %.2f Magnitude)', peak_frequency_cheby, smoothed_fft_u_cheby(peak_index_cheby)), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'r');
xlabel('log of Frequency (Hz)');
ylabel('log of Magnitude');
title('Smoothed FFT (Chebyshev)');
grid on;
hold off;

% Hann window result
subplot(2, 2, 4);
loglog(frequencies, abs(smoothed_fft_u_hann(1:length(frequencies))), 'm-', 'LineWidth', 1.5);
hold on;
loglog(frequencies, frequencies.^(-5/3), 'k--', 'LineWidth', 1.5);
scatter(peak_frequency_hann, abs(smoothed_fft_u_hann(peak_index_hann)), 'r', 'filled');
text(peak_frequency_hann, abs(smoothed_fft_u_hann(peak_index_hann)), ...
    sprintf('(%.2f Hz, %.2f Magnitude)', peak_frequency_hann, abs(smoothed_fft_u_hann(peak_index_hann))), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'r');
xlabel('log of Frequency (Hz)');
ylabel('log of Magnitude');
title('Smoothed FFT (Hann Window)');
grid on;
hold off;

sgtitle('Comparison of Different Smoothing Techniques for FFT of u - Log Plot');

%%
% 98.9403*a^7 - 402.8203*a^6 + 672.0769*a^5 - 592.7261*a^4 + 297.2906*a^3 - 87.1252*a^2 + ...
% 13.5503*a - 0.0144;
syms a;
a_original = a^9 + a*0.1 ;

a_analytical = matlabFunction(a_original,'Vars',a);

% define 3 1D grid spacings
a_spacing = linspace(0,1.1,100);
a_spacing2 = linspace(0,1.1,50);
a_spacing3 = linspace(0,1.1,10);

%define 3 values based on spacing
a_test_values = a_analytical(a_spacing);
a_test_values2 = a_analytical(a_spacing2);
a_test_values3 = a_analytical(a_spacing3);

figure(4);
plot(a_spacing, subs(a_original,a,a_spacing), 'go', 'MarkerSize',2);
xlabel('$x/\delta$','Interpreter','latex');
ylabel('w - mean velocity','Interpreter','latex')
hold on;
plot(a_spacing,a_test_values,'--k');
plot(a_spacing2,a_test_values2,'--r');
plot(a_spacing3,a_test_values3,'--b');
xlim([0 1.1]);
ylim([0 1.2])
title('Grid Testing')
legend('original','analytical - 100','analytical - 50','analytical - 10',...
    'Location','northwest')
grid on;