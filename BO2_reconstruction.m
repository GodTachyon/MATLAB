clear all;
close all;
clc;

%% generating sampled signal
nt = 51;
t = linspace(0,2*pi,nt);
f = sin(t);

%% reconstruction
% Fit a polynomial of the specified degree to the data
degree = 9;
coefficients = polyfit(t, f, degree);

% Create a range of x values for the fitted polynomial curve
nt_recon = 1001;
t_fit = linspace(min(t), max(t), nt_recon);

% Evaluate the polynomial at the x_fit values
f_fit = polyval(coefficients, t_fit);


%% visualisation
nt = 51;
t2 = linspace(0,4*pi,nt);
f2 = sin(t2);

nt_recon = 1001;
t_fit2 = linspace(min(t2), max(t2), nt_recon);

% Evaluate the polynomial at the x_fit values
f_fit2 = polyval(coefficients, t_fit2);

figure(1)
plot(t,f,'ko','MarkerSize',10,'MarkerFaceColor','k')
hold on;
plot(t_fit2,f_fit2,'--k','LineWidth',2)
plot(t2,f2,'-r','LineWidth',2)
xlim([0,max(t2)])
ylim([-1.1,1.1])
