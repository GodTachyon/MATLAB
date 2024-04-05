clear all; close all; clc;

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
%ST_ave_dat = xlsread('Reynolds_stresses.xlsx','Reynolds_stresses');

% Read data from the Excel sheet
[numData, textData, rawData] = xlsread('Reynolds_Stresses.xlsx');

% Extract data into MATLAB variables
x = numData(:, 1);
w = numData(:, 2);
ww = numData(:, 3);
uu = numData(:, 4);
vv = numData(:, 5);
uw = numData(:, 6);

figure(1);
plot(w,x);
ylim([0,1]);
ylabel('$x/\delta$','Interpreter','latex','FontSize',14);
xlabel('w velocity','Interpreter','latex','FontSize',14);




