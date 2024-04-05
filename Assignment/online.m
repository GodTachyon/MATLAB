clear all
close all

FileName = 'time_samples.hdf5';

File = h5info(FileName);

% Sampling time
% t_smpl = h5read(h5info().GroupHierarchy.Datasets.(1))
t_smpl = h5read(File.Filename,'/t');

% Wall-normal location of the samples
% x_smpl = hdf5read(hinfo().GroupHierarchy.Datasets(5)) + 1.0;
x_smpl = h5read(File.Filename,'/y');

% Sampled velocity components
% w_smpl = hdf5read(hinfo().GroupHierarchy.Datasets(2));
w_smpl = h5read(File.Filename,'/w');
% u_smpl = hdf5read(hinfo().GroupHierarchy.Datasets(3));
u_smpl = h5read(File.Filename,'/u');
% v_smpl = hdf5read(hinfo().GroupHierarchy.Datasets(4));
v_smpl = h5read(File.Filename,'/v');

% Find the index corresponding to x=0.06
x_index = find(x_smpl <= 0.06, 1, "last");

W_x_p06 = interp1(x_smpl, w_smpl.', 0.06, 'linear', 'extrap').';
W_x_n06 = interp1(x_smpl, w_smpl.', -0.06, 'linear', 'extrap').';

figure                                                          % Explore Data
surfc(x_smpl, t_smpl, w_smpl, 'EdgeColor','interp')
hold on
plot3(0.06*ones(size(t_smpl)), t_smpl, W_x_p06, '-r', 'LineWidth',1)
plot3(-0.06*ones(size(t_smpl)), t_smpl, W_x_n06, '-g', 'LineWidth',3)
hold off
grid
colormap(turbo)
xlabel('x\_smpl')
ylabel('t\_smpl')
zlabel('w\_smpl')
view(60, 30)

% t_stats = [mean(diff(t_smpl))  std(diff(t_smpl))]

Ts = mean(diff(t_smpl));
Fs = 1/Ts;
Fn = Fs/2;
L = numel(t_smpl);NFFT = 2^nextpow2(L);
FTw = fft((W_x_n06 - mean(W_x_n06)).*hann(L), NFFT)/L;
Fv = linspace(0, 1, NFFT/2+1)*Fn;
Iv = 1:numel(Fv);

[FTw_max, idx] = max(abs(FTw(Iv))*2);

figure
plot(Fv, abs(FTw(Iv))*2)
grid
xlabel('Frequency')
ylabel('Magnitude')
title('Fourier Transform')
xlim([0 6])
text(Fv(idx), FTw_max, sprintf('\\leftarrow Magnitude: %.6f\n     Frequency: %.6f',FTw_max,Fv(idx)), 'Vert','top')


% Extract wall-normal velocity signal at x=0.06
w_signal = w_smpl(:, x_index);

% Ensure the vectors have the same length
min_length = min(length(t_smpl), length(w_signal));
t_smpl = t_smpl(1:min_length);
w_signal = w_signal(1:min_length);

% Perform FFT analysis
delta_t = mean(diff(t_smpl)); % Average time step
Fs = 1 / delta_t; % Sampling frequency
L = length(w_signal);
frequencies = Fs * (0:(L-1)) / L; % Corrected the frequency calculation
fft_result = fft(w_signal);

% Ensure fft_result has valid length
valid_length = min(L, floor(L/2));
fft_amplitude = 2 / L * abs(fft_result(1:valid_length));

% Plot default FFT spectrum
figure;
subplot(2, 1, 1);
plot(t_smpl, w_signal);
title('Wall-normal Velocity Signal at x=0.06');
xlabel('Time (s)');
ylabel('Velocity (m/s)');


subplot(2, 1, 2);
plot(frequencies, fft_amplitude, 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('FFT Spectrum - Default Plot');
grid on;
hold off;


% Logarithmic plot
figure;
semilogx(frequencies, fft_amplitude, 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('FFT Spectrum - Logarithmic Plot');
grid on;
hold off;


% Apply window functions and plot FFT spectra
window_functions = {'rectwin', 'hamming', 'hann', 'blackman'};
figure;
for i = 1:length(window_functions)
window_function = window_functions{i};
window = str2func(window_function);

% Apply window function
w_windowed = w_signal .* window(L)';

% Perform FFT analysis on windowed signal
fft_result_windowed = fft(w_windowed);
fft_amplitude_windowed = 2/L * abs(fft_result_windowed(1:valid_length));

% Plot FFT spectrum with window function
subplot(length(window_functions), 1, i);
plot(frequencies, fft_amplitude_windowed, 'LineWidth', 2, 'DisplayName', window_function);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title(['FFT Spectrum with ' window_function ' Window']);
grid on;
hold off;


end

xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('FFT Spectrum with Different Window Functions');
legend('show');
grid on;
hold off;



% Moving average filtering
window_size = 10; % Adjust window size as needed
w_smoothed = movmean(w_signal, window_size);

% Plot original and smoothed signals
figure;
subplot(2, 1, 1);
plot(t_smpl, w_signal, 'LineWidth', 2, 'DisplayName', 'Original Signal');
title('Original and Smoothed Signals');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('show');
grid on;
hold off;


subplot(2, 1, 2);
plot(t_smpl, w_smoothed, 'r--', 'LineWidth', 2, 'DisplayName', 'Smoothed Signal');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('show');
grid on;
hold off;