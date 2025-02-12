% Define the mass matrix (example)
M = [0 0; 0 1]; % Mass matrix for a simple DAE system

% Define the right-hand side function
f = @(t, y) [144*y(1)-7; -y(2)]; % Example RHS function

% Define the Jacobian function (for the ODE part)
J = @(t, y) [144 0; 0 -1]; % Jacobian of f(t, y) with respect to y

% Set options with the mass matrix and Jacobian
options = odeset('Mass', M, 'Jacobian', J);

% Initial conditions
y0 = [0.49; 2];

% Time span
tspan = [0 10];

% Solve the DAE system
[t, y] = ode15s(f, tspan, y0, options);

plot(t, y)