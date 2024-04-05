%% COMPUTE STATS
% This file has been created at Cranfield University in order to
% provide supplementary materials for the Data Analysis and Uncertainty
% module assignment (N-CFD-DAPP)
%
% The script reads data obtained from direct numerical simulation of
% a low Reynolds number channel flow. It reads the excel data sheet
% named "Reynolds_stresses.xlsx" including spatio-temporal averages 
% computed from 200,000 time steps. The data sheet also contains the 
% parameters of the channel flow.
%
% Spatio-temporal averages are computed by computing the time-averaged 
% value at every cell thorughout the simulations and then computing the
% spatial averages of the time-averaged data along planes parallel to the
% wall. Spatial-averaging is denoted by <> in this script.
%
% The script also reads the 'time_samples.hdf5' file which contains time
% samples of the velocity components and the pressure along a wall-normal
% line from the bottom wall to the centre line of the channel. This data is
% processed using time-averaging symbolised by <>_t.
%
% The data is made non-dimensional based on the fluid density, the average
% streamwise velocity, and the channel half-height.
%
% For equations Kundu's book has been used (fifth edition)
%
% written by Tamás István Józsa
% 21/09/2023
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