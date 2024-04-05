% Code for calculating the wall spacing based on y+ value. 
% Equations taken from Cadence online calculator.
% Developed by Aswath Ashok on 16/03/2024 at Cranfield University.

% Input values (All in SI units)
Re = 40000; % Reynolds number (-)
y_plus = 30; % Desired y+ value (-)
dyn_viscosity = 7.87e-6; % Dynamic viscosity (kg/m.s)
len = 0.07; % Characteristic length scale (m)
density = 1; % Density (kg/m^3)

% Calculation
velocity = (Re*dyn_viscosity)/(density*len); % velocity of flow (m/s)
c_f = 0.026/(Re^(1/7)); % Friction coeffcient (-)
tau_w = (c_f*density*velocity^2)/2; % Wall friction factor (kg/m.s^2)
velocity_fric = sqrt(tau_w/density); % Friction velocity (m/s)
wall_space = (y_plus*dyn_viscosity)/(velocity_fric*density); % Wall spacing for mesh (m)

