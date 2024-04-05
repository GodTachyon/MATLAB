figure;

% Histogram and Fitted Distribution subplot
subplot(1, 2, 1);
hhf = histfit(u_subset, 2.5E+2, 'normal');
xd  = hhf(2).XData;
yd  = hhf(2).YData;
parms = fitdist(u_subset,'normal');
xv = [-1 0 1]*1.96*parms.sigma+parms.mu;
yv = interp1(xd, yd, xv);
hold on
hp1 = plot([1;1]*xv(2), [0;1]*yv(2), '--r', 'LineWidth',2.5, 'DisplayName','Mean');
hp2 = plot([1;1]*xv([1 3]), [0;1]*yv([1 3]), ':r', 'LineWidth',2.5, 'DisplayName','95% Limits');
hold off
grid
xlabel('‘u\_subset’ Amplitude')
ylabel('Frequency')
legend([hp1; hp2(1)], 'Location','best')
title('Histogram and Fitted Distribution');

% PDF subplot
subplot(1, 2, 2);
x_values = linspace(min(u_subset), max(u_subset, 100));
pdf_values = pdf(parms, x_values);
plot(x_values, pdf_values, 'b', 'LineWidth', 2);
hold on;

% Marking the mean on the PDF
mean_value = mean(u_smpl);
interp_mean = interp1(x_values, pdf_values, mean_value, 'linear', 'extrap');
plot(mean_value, interp_mean, 'ro', 'MarkerSize', 8, 'DisplayName', 'Mean');

hold off;
grid;
xlabel('‘u\_subset’ Amplitude');
ylabel('Probability Density');
legend('PDF', 'Mean', 'Location', 'best');
title('PDF of u\_subset with Mean');

%%
% Assuming you have a function f representing the differential equation
f = @(x) x.^2 - 1;

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
    numerical_values = zeros(size(x_values));  % Replace with actual computation
    
    % Compute error
    errors(i) = norm(numerical_values - a_analytical(x_values), inf);
end
% Perform Richardson extrapolation - log
p = polyfit(log(1./grid_sizes), log(errors), 1);
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


