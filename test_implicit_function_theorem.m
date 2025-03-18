syms x1 x2 y1 y2
F1 = x1 + 1*y1*y2 + 1*y2;  % First equation
F2 = x2 + y1 + 1*y2; % Second equation
F3 = x1 - y2 - cos(y1);
F4 = x1 + sin(y1) - y2;

% Define F as a vector
F = [F3; F4];

% Define variables
X = [x1; x2];  % Independent variables
Y = [y1; y2];  % Dependent variables

% Compute Jacobian of F with respect to Y
J = jacobian(F, Y);

% Check if the Jacobian is invertible
det_J = det(J);

disp('Jacobian matrix:')
disp(J)
disp('Determinant of Jacobian:')
disp(det_J)

% If det(J) ≠ 0, solve for Y in terms of X
if det_J ~= 0
    G = solve(F, Y); % Solve for y in terms of x
    g_x = [G.y1(1); G.y2(1)];
    disp('Inverse function g(x):')
    disp(g_x)
else
    disp('Jacobian is singular; Implicit Function Theorem does not guarantee an inverse.')
end
