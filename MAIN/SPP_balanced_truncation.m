addpath(genpath('../../jti-code'))

load("modelSPP_no_crosstalk_full.mat")

I = model.I;
par = model.par;
X0 = model.X0;

A0 = X0(I.A);
S0 = X0(I.A);

% B C D
A = [-par(I.k_b), 0, 0;
     0, -par(I.k_c), 0;
     par(I.k_bd), par(I.k_cd), -par(I.k_d)];

B = [S0*par(I.k_ab);
     S0*par(I.k_ac);
     0];

C = [0, 0, 1];

D = 0;

% Create state-space model
model_ss = ss(A, B, C, D);
fprintf('Original system order: %d\n', order(model_ss));

% Check poles and zeros
poles = pole(model_ss);
zeros = zero(model_ss);

fprintf('Poles: '); disp(poles');
fprintf('Zeros: '); disp(zeros');

% Cancel poles and zeros
model_ss_minimal = minreal(model_ss, 1e-6);

% Compute Hankel singular values to see the "importance" of each state
hsv = hsvd(model_ss_minimal);
fprintf('Hankel singular values: ');
disp(hsv');

% Apply balanced truncation
reduced_order = 2;
model_ss_red = balred(model_ss_minimal, reduced_order);

% plot
f = figure;

t = 0:0.001:0.5;
u = 25 * ones(size(t));

hold on
lsim(model_ss, u, t);
lsim(model_ss_red, u, t);

hold on
[y, t] = step(model_ss);
semilogy(t, y);
hold on
[y, t] = step(model_ss_red);
semilogy(t, y);
hold on