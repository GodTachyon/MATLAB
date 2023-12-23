% The FFT.m file reads the DNS data and now manipulate the wall normal
% velocity component (u) at X = 0.06 location. FFT is used to analyse the
% magnitude - frequency relationship. Effects of regular plot and log-log
% plot is also investigated. Filtering methods are utilised and their
% effects explored.
% Statistical methods are then carried out to find out the effectiveness of
% wall normal velocity component (u). Mean, Standard Deviation, Skewness
% and Kurtosis quantities are utilised in creating a Probability Density
% Function (PDF) and plotted. Pearson Correlation is utilised to test
% relationship between the u - component and v - component.
%
% Code for Reading the initial data was provided by Dr. Tamas Jozsa for the
% Data and Uncertainity Analysis Module.
% Further development was done by Aswath Ashok for the Assignment for the
% above module.
% Created on 22/12/23. Updates can be checked on github.

%% READ CHANNEL PARAMETERS AND SPATIO-TEMPORAL AVERAGES

% read parameters
% Lz, Lx, Ly, nu, Delta p
params = xlsread('Reynolds_stresses.xlsx','parameters');

Lz = params(1); Lx = params(2); Ly = params(3);
nu = params(4); % kinematic viscosity

% these values are equal to unity because they are the reference quantities
% used to make the data dimensionless
u_b = 1.0; % bulk velocity (average velocity in the entire channel)
rho = 1.0; % density
delta = Lx/2; % boundary layer thickness  =  channel half-height

% bulk Reynolds number based on channel half height and mean velocity
Re_b  =  u_b*delta/nu;

% read wall-normal coordinate and spatio-temporal averages
% x, <w>, <w'w'>, <u'u'> , <v'v'>, <u'w'>
ST_ave_dat = xlsread('Reynolds_stresses.xlsx','Reynolds_stresses');


%% READ TIME SAMPLES AT PROBES PLACED ALONG A WALL_NORMAL LINE

hinfo  =  hdf5info('time_samples.hdf5');

% sampling time
t_smpl = hdf5read(hinfo.GroupHierarchy.Datasets(1));
% wall-normal location of the samples
x_smpl = hdf5read(hinfo.GroupHierarchy.Datasets(5))+1.0;

% sampled velocity components
% each row represents a time instant as dictated by t_smpl
% each column represents a spatial location as dictated by y_smpl
w_smpl = hdf5read(hinfo.GroupHierarchy.Datasets(2));
u_smpl = hdf5read(hinfo.GroupHierarchy.Datasets(3));
v_smpl = hdf5read(hinfo.GroupHierarchy.Datasets(4));


%% instantaneous velocity plots
figure(1)
hold on
plot(x_smpl,u_smpl(1,:),':r',LineWidth=2)
plot(x_smpl,v_smpl(1,:),'--b',LineWidth=2)
%plot(x_smpl,w_smpl(1,:),'-k',LineWidth=2)
xlabel('$x/\delta$','Interpreter','latex','FontSize',14)
xlim([0,1])
legend('$u/u_b$','$v/u_b$','$w/u_b$',...
    'Interpreter','latex','FontSize',14,'Location','east')
set(gca,'TickLabelInterpreter','latex','FontSize',14)
hold off

%% FFT TIME BOIIIIII

% Extract a subset of data
u_subset = u_smpl(1:200000,34);
v_subset = v_smpl(1:200000,34);

% Calculate the FFT of the functions
fft_u = fft(u_subset);
fft_v = fft(v_subset);

% Frequency vector (assuming the time steps are in seconds)
fs = 1/(t_smpl(2)-t_smpl(1));  % Sampling frequency
frequencies = linspace(0, fs/2, length(fft_u)/2 + 1);

% Moving Average filter
windowSize = 2;  % Adjust as needed
smoothed_fft_u = smoothdata(abs(fft_u), 'movmean', windowSize);

% Plot the magnitude spectrum of u
figure(2);
subplot(2,1,1);
plot(frequencies,abs(fft_u(1:length(frequencies))));
xlim([0,5]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT of u');
% Plot the magnitude spectrum of v
% subplot(2,1,2);
% plot(frequencies,abs(fft_v(1:length(frequencies))));
% xlim([0,10]);
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% title('FFT of v');

% Hamming Filter (MAKE IT WORK SOMEHOW!!!)
subplot(2,1,2);
plot(frequencies, smoothed_fft_u(1:length(frequencies)));
xlim([0,5]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Windowed FFT (moving average)');

% Plot both the original and smoothed FFT on a log-log scale
figure (3);
loglog(frequencies, abs(fft_u(1:length(frequencies))), 'b-', 'LineWidth', 1.5);
hold on;
loglog(frequencies, smoothed_fft_u(1:length(frequencies)), 'r-', 'LineWidth', 1.5);
% Plot a line with a slope of -5/3
slope_line = frequencies.^(-5/3);
loglog(frequencies, slope_line, 'g--', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ylim([10^-3,10^3]);
title('FFT of u and Smoothed FFT (Moving Average) - Log-Log Scale');
legend('Original FFT', 'Smoothed FFT','Energy Cutoff');
hold off;

%% Probability Density Function TIMEEE BOIIIIIII

% Calculate the PDF of u using the pdf function
numPoints = 126;  % Number of points for PDF estimation
n_bins = 126; % Number of points for Histogram
x = linspace(-0.3, 0.2, numPoints);
u_mean = mean(u_smpl);
u_std = std(u_smpl);
u_skew = skewness(u_smpl);
u_kurtosis = kurtosis(u_smpl);
pdfValues = pdf('Normal', x, u_mean,u_std, u_skew, u_kurtosis);

% Plot the PDF of u
figure (4);
subplot(1,2,1)
plot(x, pdfValues, 'b-', 'LineWidth', 1.5);
set(gca,'XLim',[-0.3 0.3],'XTick',-0.3:0.10:0.3)
xlabel('Velocity (u)');
ylabel('Probability Density');
title('Probability Density Function (PDF) of Velocity u');
% Histogram Plot of wall normal velocity i.e. u
subplot(1,2,2)
histogram(u_smpl, 'Normalization','pdf');
set(gca,'XLim',[-0.2 0.2],'XTick',-0.2:0.10:0.2)
xlabel('Velocity (u)');
ylabel('Probability Density');
title('Histogram of the PDF of Velocity u');

% plot skew and kurtosis
% figure (5);
% hold on
% plot(x, u_skew);
% plot(x,u_kurtosis);

%% Pearson Correlation

% Calculate the Pearson correlation coefficient
correlation_coefficient = corrcoef(u_smpl, v_smpl);
pearson_correlation = correlation_coefficient(1, 2);

% Plot the data
figure (6);
hold on;
scatter(u_smpl, v_smpl, 'b.');
% Plot y = x line
plot(get(gca, 'XLim'), get(gca, 'XLim'), 'r--', 'LineWidth', 1.5);
xlabel('X Velocity (u)');
ylabel('Y Velocity (v)');
title(['Scatter Plot of u and v - Pearson Correlation: ' num2str(pearson_correlation)]);

% Display the Pearson correlation coefficient
disp(['Pearson Correlation Coefficient: ' num2str(pearson_correlation)]);





