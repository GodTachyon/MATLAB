% Test for mean plots.
% 
% Code written by Aswath Ashok for the Data Analysis module
% at Cranfield University for the MSc CFD course.

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
% figure(1);
% yyaxis left
% plot(x_smpl,w_smpl(1,:),'-k',LineWidth=2);
% ylabel('Normalized Velocity of w','Interpreter','latex','FontSize',14);
% yyaxis right
% plot(x_smpl,u_smpl(1,:),'-r',LineWidth=2);
% hold on;
% plot(x_smpl,v_smpl(1,:),'-b',LineWidth=2);
% xlabel('$x/\delta$','Interpreter','latex','FontSize',14);
% ylabel('Normalized Velocity of u and w','Interpreter','latex','FontSize',14);
% xlim([0,1]);
% legend('$w/u_b$','$u/u_b$','$v/u_b$',...
%     'Interpreter','latex','FontSize',14,'Location','east');
% set(gca,'TickLabelInterpreter','latex','FontSize',14);
% hold off;
% Upto here is from the provided assignment file.



%% MEAN OF VELOCITY COMPONENTS FOR 126 Sampling points
u_mean = mean(u_smpl);
v_mean = mean(v_smpl);
w_mean = mean(w_smpl);

figure(2);
yyaxis left
plot(x_smpl,w_mean,'-k',LineWidth=2);
ylabel('Mean Velocity w','Interpreter','latex','FontSize',14);
yyaxis right
plot(x_smpl,u_mean,'-r',LineWidth=2);
hold on
plot(x_smpl,v_mean,'-b',LineWidth=2);
ylabel('Mean Velocity u and v','Interpreter','latex','FontSize',14);

xlabel('$x/\delta$','Interpreter','latex','FontSize',14);
legend('w-mean','u-mean','v-mean','FontSize',14,'Location','southeast');
title('Mean of Velocity components','Interpreter','latex','FontSize',14);
set(gca,'TickLabelInterpreter','latex','FontSize',10);
grid on
hold off;


%% STANDARD DEVIATION OF THE VELOCITY COMPONENTS

u_std = std(u_smpl);
v_std = std(v_smpl);
w_std = std(w_smpl);

figure(3);
plot(x_smpl,u_std,'-r',LineWidth=2);
hold on;
plot(x_smpl,v_std,'-b',LineWidth=2);
plot(x_smpl,w_std,'-k',LineWidth=2);
xlabel('$x/\delta$','Interpreter','latex','FontSize',14);
ylabel('STD Velocity','Interpreter','latex','FontSize',14);
xlim([0,1]);
title('Standard Deviation of Velocity Components','Interpreter','latex','FontSize',14);
legend('u-std','v-std','w-std',...
    'Interpreter','latex','FontSize',14,'Location','east');
set(gca,'TickLabelInterpreter','latex','FontSize',10);
grid on
hold off;

%% Reynolds_stresses.xlsx spatio-temporal values plot

y = ST_ave_dat(:,1); w = ST_ave_dat(:,2); ww = ST_ave_dat(:,3); 
uu = ST_ave_dat(:,4); vv = ST_ave_dat(:,5); uw = ST_ave_dat(:,6);

% figure(4);
% subplot(1,2,1);
% plot(y,w,'--b',LineWidth=2);
% xlabel('X','FontSize',14);
% ylabel('spatio-temporal averages','FontSize',14);
% xlim([0,1]);
% legend('w-avg');
% subplot(1,2,2);
% plot(y,ww,':r',LineWidth=2);
% hold on;
% plot(y,uu,'--b',LineWidth=2);
% plot(y,vv,'-k',LineWidth=2);
% plot(y,uw,'-g',LineWidth=2);
% xlabel('X','FontSize',14);
% ylabel('spatio-temporal averages','FontSize',14);
% xlim([0,1]);
% legend('ww-avg','uu-avg','vv-avg','uw-avg');
% hold off;

%% Gradient of u_mean, v_mean and w_mean

% Spacing for each componenet

velocity_spacing = x_smpl(2) - x_smpl(1);

u_grad = gradient(u_mean,velocity_spacing);
v_grad = gradient(v_mean,velocity_spacing);
w_grad = gradient(w_mean,velocity_spacing);


figure(5);
% Plot w on the left y-axis
yyaxis left;
plot(x_smpl, w_grad, '-k', 'LineWidth', 1.5);
ylabel('w Component');
% Plot u and v on the right y-axis
yyaxis right;
plot(x_smpl, u_grad, '-r', 'LineWidth', 1.5);
hold on;
plot(x_smpl, v_grad, '-b', 'LineWidth', 1.5);
ylabel('u and v Components');

% Draw a vertical line at x = 0.06 for all components
%xline(0.06, 'k--', 'LineWidth', 1.5);

% Add labels and title
xlabel('$x/\delta$','Interpreter','latex','FontSize',14);
title('Gradient of Velocity Components','Interpreter','latex','FontSize',14);
% Add a legend
legend('w - grad', 'u - grad', 'v - grad','Interpreter','latex','FontSize',14);
set(gca,'TickLabelInterpreter','latex','FontSize',10);
grid on
% Add text label for x = 0.06 on the x-axis
% text(0.06, min(ylim), 'x = 0.06', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
%     'FontSize', 10, 'FontWeight', 'bold')

