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

%% MEAN OF VELOCITY COMPONENTS FOR 126 Sampling points
u_mean = mean(u_smpl);
v_mean = mean(v_smpl);
w_mean = mean(w_smpl);

% figure(1);
% yyaxis left
% plot(x_smpl,w_mean,'-k',LineWidth=2);
% ylabel('Mean Velocity w','Interpreter','latex','FontSize',14);
% yyaxis right
% plot(x_smpl,u_mean,'-r',LineWidth=2);
% hold on
% plot(x_smpl,v_mean,'-b',LineWidth=2);
% ylabel('Mean Velocity u and v','Interpreter','latex','FontSize',14);
% 
% xlabel('$x/\delta$','Interpreter','latex','FontSize',14);
% legend('w-mean','u-mean','v-mean','FontSize',14,'Location','southeast');
% title('Mean of Velocity components');
% set(gca,'TickLabelInterpreter','latex','FontSize',10);
% grid on
% hold off;

%% Gradient of u_mean, v_mean and w_mean

% Spacing for each componenet
velocity_spacing = x_smpl(2) - x_smpl(1);

u_grad = gradient(u_mean,velocity_spacing);
v_grad = gradient(v_mean,velocity_spacing);
w_grad = gradient(w_mean,velocity_spacing);


% figure(2);
% % Plot w on the left y-axis
% yyaxis left;
% plot(x_smpl, w_grad, '-k', 'LineWidth', 1.5);
% ylabel('w Component');
% % Plot u and v on the right y-axis
% yyaxis right;
% plot(x_smpl, u_grad, '-r', 'LineWidth', 1.5);
% hold on;
% plot(x_smpl, v_grad, '-b', 'LineWidth', 1.5);
% ylabel('u and v Components');
% 
% % Add labels and title
% xlabel('X');
% title('Gradient of Velocity Components');
% % Add a legend
% legend('w - grad', 'u - grad', 'v - grad');


%% Use symbolic to interpolate??

syms a;

% Define the original function
a_original = a^9 + a*0.1 ;

% Convert the symbolic function to a MATLAB function
a_analytical = matlabFunction(a_original,'Vars',a);

% Reference value at x = 0.2
reference_value = a_analytical(0.2);

% Define grid spacings
grid_sizes = [100, 50, 10];
errors = zeros(size(grid_sizes));

% Perform Richardson extrapolation and calculate the grid convergence index
refinement_factors = zeros(size(grid_sizes));
grid_convergence_indices = zeros(size(grid_sizes));

for i = 1:length(grid_sizes)
    % Define grid spacing
    a_spacing = linspace(0, 1.1, grid_sizes(i));
    
    % Calculate values based on spacing
    a_test_values = a_analytical(a_spacing);
    
    % Compute error using Richardson extrapolation
    h = a_spacing(2) - a_spacing(1);
    errors(i) = richardsonExtrapolation(reference_value, a_test_values(2), h);
    
    % Calculate refinement factor and grid convergence index
    if i > 1
        refinement_factors(i) = errors(i-1) / errors(i);
        grid_convergence_indices(i) = log2(refinement_factors(i));
    end
    
    % Plot the results
    figure;
    plot(a_spacing, subs(a_original, a, a_spacing), 'go', 'MarkerSize', 2);
    hold on;
    plot(a_spacing, a_test_values, '--k');
    xlim([0 1.1]);
    ylim([0 1.2]);
    legend('original', ['analytical - ' num2str(grid_sizes(i))]);
    title(['Grid Convergence Study - Grid Size: ' num2str(grid_sizes(i))]);
    xlabel('a');
    ylabel('Function Value');
end

% Calculate asymptotic range of GCI values
asymptotic_range = grid_convergence_indices(end) - grid_convergence_indices;

% Display a table of grid sizing, refinement, error, and asymptotic range
table_data = table(grid_sizes', refinement_factors', errors', grid_convergence_indices', asymptotic_range', 'VariableNames', {'Grid Size', 'Refinement Factor', 'Error', 'Grid Convergence Index', 'Asymptotic Range'});
disp('Grid Convergence Study Results:');
disp(table_data);




%% Richardson Time BOIIIII

num_refinements = 10;
%  Initialize arrays for errors and grid sizes
errors = zeros(num_refinements, 1);
grid_sizes = zeros(num_refinements, 1);

%  Perform grid convergence analysis
for i = 1:num_refinements
    % Refine the grid
    grid_size =  2*i;
    grid_sizes(i) = grid_size;
    x_values = linspace(0, 1.1, grid_size);
    
    % Numerical solution (replace with your numerical method)
    numerical_values = abs(zeros(x_values));  % Replace with actual computation
    
    % Compute error
    errors(i) = norm(numerical_values - a_analytical(x_values), inf);
end
% Perform Richardson extrapolation - log
p = polyfit(log(1./grid_sizes), log(errors), 2);
converged_error = exp(polyval(p, log(0)));
% richardson normal
p_normal = polyfit((1./grid_sizes), errors,1);
converged_normal = polyval(p_normal,0);

figure (5);
plot(1./grid_sizes, errors, '-o');
hold on;
plot(1./grid_sizes, converged_normal, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
%  Verify convergence by plotting error vs grid size
figure (6);
loglog(1./grid_sizes, errors, '-o');
hold on;
loglog(0, converged_error, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Grid Size (1/dx)');
ylabel('Error');
title('Grid Convergence Analysis');
legend('Error vs Grid Size', 'Converged Error');
grid on;
hold off;


function error = richardsonExtrapolation(reference_value, refined_value, h)
    % Richardson Extrapolation formula
    error = abs(reference_value - refined_value) / (2^2 - 1);
end
