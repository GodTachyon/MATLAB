%% Polyfit to find the above symbolic function

% Given data points
x = x_smpl;
y = w_mean;

% Degree of the polynomial (you can adjust this based on the desired degree)
degree = 7;

% Find the coefficients of the best-fitting polynomial
coefficients = polyfit(x, y, degree);


% Generate a polynomial using the coefficients
bestFitPoly = polyval(coefficients, x);

% Plot the original data points and the best-fitting polynomial
figure(1);
plot(x, y, 'o', x, bestFitPoly, '-');
xlabel('x');
ylabel('y');
title('Best-Fitting Polynomial');
legend('Data Points', 'Best-Fitting Polynomial');

% Display the coefficients of the polynomial
% disp('Coefficients of the best-fitting polynomial:');
% disp(coefficients);

%% try different spacing

% Symbolic function 9th degree
syms a b;
% a_original = 30.4552*a^5 - 87.1451*a^4 + 92.7845*a^3 - 45.3534*a^2 + ...
%     10.6033*a + 0.0298;
% a_analytical = matlabFunction(a_original,'Vars',a);
% a_spacing = linspace(0,1.1,50);
% a_test_values = a_analytical(a_spacing);
% 
% figure(2);
% plot(a_spacing, subs(a_original,a,a_spacing), 'ko', 'MarkerSize',2);
% ylabel('original - mean')
% hold on;
% plot(a_spacing,a_test_values,'--b');
% xlim([0 1.5]);

b_test = sin(1/3*pi*b);
b_diff = diff(b_test,b);
b_analytical = matlabFunction(b_diff,'Vars',b);
b_spacing = linspace(0,1/3*pi,10);
b_test_values = b_analytical(b_spacing);

figure(2);
plot(b_spacing, subs(b_test,b,b_spacing), 'ko');
hold on;
plot(b_spacing,b_test_values,'--b');
xlim([0 1.5]);



